import pandas as pd
import numpy as np
import re


# Data preprocessing
## GRS1
### Helper functions
def classify_genotype(row):
    DR3= 'DRB1*03:01-DQA1*05:01-DQB1*02:01'
    DR4DQ8= 'DRB1*04:01-DQA1*03:01-DQB1*03:02'

    h1 = str(row["haplotype1"])
    h2 = str(row["haplotype2"])

    is_dr3_1 = h1 == DR3
    is_dr3_2 = h2 == DR3

    is_dr4_1 = h1 == DR4DQ8
    is_dr4_2 = h2 == DR4DQ8

    # Both haplotypes are DR3
    if is_dr3_1 and is_dr3_2:
        return "DR3/DR3"

    # Both haplotypes are DR4DQ8
    elif is_dr4_1 and is_dr4_2:
        return "DR4-DQ8/DR4-DQ8"

    # One DR3 and one DR4DQ8
    elif (is_dr3_1 and is_dr4_2) or (is_dr4_1 and is_dr3_2):
        return "DR3/DR4-DQ8"

    # One haplotype is DR3
    elif is_dr3_1 or is_dr3_2:
        return "DR3/X"

    # One haplotype is DR4DQ8
    elif is_dr4_1 or is_dr4_2:
        return "DR4-DQ8/X"

    # Neither haplotype matches DR3 or DR4DQ8
    else:
        return None

### Data preparation
train_hla_df = pd.read_csv('./GRS-input/Discovery_HLA-haplotype.csv').set_index('ID')

train_hla_df["genotype"] = train_hla_df.apply(classify_genotype, axis=1)

train_wide_df = (
    train_hla_df
    .pivot_table(
        index="genotype",
        columns=train_hla_df.index,
        aggfunc="size",
        fill_value=0
    )
    .reindex(columns=train_hla_df.index, fill_value=0)
)

val_hla_df = pd.read_csv('./GRS-input/Validation_HLA-haplotype.csv').set_index('ID')

val_hla_df["genotype"] = val_hla_df.apply(classify_genotype, axis=1)

val_wide_df = (
    val_hla_df
    .pivot_table(
        index="genotype",
        columns=val_hla_df.index,
        aggfunc="size",
        fill_value=0
    )
    .reindex(columns=val_hla_df.index, fill_value=0)
)

## GRS2 and CGRS

### Non-DR-DQ loci that can use genotype based data in GRS2

beta_df = pd.read_csv("../GRS/GRS2/DC181785_Table_S4.csv")
beta_df = beta_df[0:16]
loci = beta_df["Locus"].dropna().unique()

def grs2_hla_matrix(df, loci):

    result = pd.DataFrame(
        0,
        index=loci,
        columns=df.index
    )

    for patient, row in df.iterrows():
        for locus in loci:
            result.loc[locus, patient] = sum(
                str(locus) in str(row[col])
                for col in df.columns
                if pd.notna(row[col])
            )

    return result

train_allele_df = pd.read_csv("./grs_input/Discovery_HLA-allele.csv",index_col=0)
val_allele_df = pd.read_csv("./grs_input/Validation_HLA-allele.csv",index_col=0)

train_allele_matrix = grs2_hla_matrix(train_allele_df, loci)
val_allele_matrix = grs2_hla_matrix(val_allele_df, loci)

## ## Construct the full feature matrix
train_matrix=pd.read_csv('./grs_new/all_snp_sample_163.csv',index_col=0)
val_matrix=pd.read_csv('./grs_new/all_snp_sample_527.csv',index_col=0)

train_matrix2=pd.read_csv('./grs_new/plus_snp_sample_163.csv',index_col=0)
val_matrix2=pd.read_csv('./grs_new/plus_snp_sample_527.csv',index_col=0)

def combine_df(a,b,c,d):
      cols = a.columns
      full_df = pd.concat([a[cols], b[cols], c[cols],d[cols]], axis=0)
      return full_df

a = pd.read_csv('./grs_new/all_snp_sample_163.csv',index_col=0)
train_matrix = combine_df(a, train_wide_df, train_allele_matrix, train_matrix2)

a= pd.read_csv('./grs_new/all_snp_sample_527.csv',index_col=0)
val_matrix = combine_df(a, val_wide_df, val_allele_matrix, val_matrix2)


## Stepwise GRS calculation

def match_pattern(pattern, value):
    """
    Check whether a beta pattern matches the patient's HLA value.

    X represents any digit
    """
    if pd.isna(pattern) or pd.isna(value):
        return False

    pattern = str(pattern)
    value = str(value)

    # Replace X with a wildcard for any digit
    regex = re.escape(pattern).replace(r'\X', r'\d')

    # A partial match is sufficient within the patient cell
    return re.search(regex, value) is not None

def step1(row, beta_df):

    h1 = str(row["haplotype1"])
    h2 = str(row["haplotype2"])

    for _, beta_row in beta_df.iterrows():

        b1 = beta_row["haplotype1"]
        b2 = beta_row["haplotype2"]

        # Forward match
        forward = (
            match_pattern(b1, h1)
            and
            match_pattern(b2, h2)
        )

        # Reverse match after swapping the two haplotypes
        reverse = (
            match_pattern(b1, h2)
            and
            match_pattern(b2, h1)
        )

        if forward or reverse:
            return beta_row["beta"]

    return 0


def match_haplotype(pattern, value):
    """Check whether the patient's HLA string contains the beta haplotype.
    X represents any digit.
    """
    if pd.isna(pattern) or pd.isna(value):
        return False

    pattern = str(pattern)
    value = str(value)

    # X -> any digit
    regex = re.escape(pattern).replace(r'\X', r'\d')

    return re.search(regex, value) is not None


def calculate_patient_beta(row,beta_df):

    total = 0.0

    h1 = str(row["haplotype1"])
    h2 = str(row["haplotype2"])

    for _, beta_row in beta_df.iterrows():

        pattern = beta_row["haplotype"]
        beta = float(beta_row["Beta"])

        # Evaluate both patient haplotypes separately
        count = 0

        if match_haplotype(pattern, h1):
            count += 1

        if match_haplotype(pattern, h2):
            count += 1

        # Add 0, 1, or 2 times the beta value depending on the count
        total += count * beta

    return total

def step1_step2(hla_df):
    """Combine the two HLA terms: use step1 first, and fall back to step2 when no interaction matches."""
    beta_df1 = pd.read_csv("../GRS/step1.csv")
    beta_df1.columns = beta_df1.columns.str.strip()
    beta_df1 = beta_df1.rename(columns={"Beta": "beta"})

    hla_df = hla_df.copy()
    hla_df["beta"] = hla_df.apply(lambda row: step1(row, beta_df1), axis=1)
    vector1 = hla_df["beta"]

    beta_df2 = pd.read_csv("../GRS/step2.csv")
    vector2 = hla_df.apply(
        lambda row: calculate_patient_beta(row, beta_df2),
        axis=1
    )

    vector = vector1.where(vector1 != 0, vector2)
    vector.name = "HLA_beta"
    return vector

## Step 3: additive SNP calculation
def snp_match(index, snp_rule):
    """
    Check whether a SNP index in the matrix matches one SNP rule in beta_df.

    . denotes any nucleotide;
    , denotes OR.
    """

    for snp in str(snp_rule).split(","):
        snp = snp.strip()

        if not snp:
            continue

        # . indicates that the last base can be any nucleotide
        if snp.endswith("."):
            pattern = re.escape(snp[:-1]) + r"[ACGT]"
        else:
            pattern = re.escape(snp)

        if re.fullmatch(pattern, str(index)):
            return True

    return False


def calculate_grs(matrix, model):

    beta_df = pd.read_csv("../GRS/step3.csv", index_col=0)

    grs = beta_df[beta_df["model"] == model].copy()

    matrix = matrix[
        ~matrix.index.duplicated(keep="first")
    ]

    grs_result = np.zeros(len(matrix.columns))

    matched_count = 0

    for snp_rule, row in grs.iterrows():

        beta = row["beta"]

        # Identify matrix rows that match the SNP rule
        matched = matrix.index[
            matrix.index.to_series().apply(
                lambda x: snp_match(x, snp_rule)
            )
        ]

        if len(matched) == 0:
            continue

        matched_count += 1

        X = matrix.loc[matched].sum(axis=0)

        grs_result += X.to_numpy() * beta

    print(
        f"Model: {model}, "
        f"matched SNPs: {matched_count}/{len(grs)}"
    )

    return pd.Series(grs_result, index=matrix.columns, name=model)


# Compute GRS scores
grs_train_all = pd.DataFrame(index=train_matrix.columns)
grs_val_all = pd.DataFrame(index=val_matrix.columns)
for model in ['GRS1', 'GRS2', 'CGRS']:
    if model == 'GRS1':
        grs_train = calculate_grs(train_matrix, model)
        grs_val = calculate_grs(val_matrix, model)
    else:
        train_vector = step1_step2(train_hla_df).reindex(train_matrix.columns)
        val_vector = step1_step2(val_hla_df).reindex(val_matrix.columns)
        grs_train = calculate_grs(train_matrix, model).add(train_vector, fill_value=0)
        grs_val = calculate_grs(val_matrix, model).add(val_vector, fill_value=0)
    grs_train_all[model] = grs_train.reindex(grs_train_all.index)
    grs_val_all[model] = grs_val.reindex(grs_val_all.index)


# Merge disease status

train_cli=pd.read_csv('./grs_input/Dis_cli.csv', index_col=0)
val_cli=pd.read_csv('./grs_input/Val_cli.csv', index_col=0)

match_df = pd.read_csv('./grs_input/Discovery_match.csv',index_col=0)

snp_to_rna = match_df.set_index("SNP_ID")["RNA_ID"]

grs_train_all["Disease"] = (
    snp_to_rna.reindex(grs_train_all.index)
    .map(train_cli["Disease"])
)

snp_to_rna = match_df.set_index("SNP_ID")["RNA_ID"]
grs_val_all['Disease'] = val_cli.loc[grs_val_all.index]['Disease']

grs_train_all.to_csv('./grs_output/grs_train_all.csv')
grs_val_all.to_csv('./grs_output/grs_val_all.csv')


