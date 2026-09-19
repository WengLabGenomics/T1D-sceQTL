library(MatrixEQTL)
arg=commandArgs(T)
ct = arg[1]

base.dir = './input/disease_all_geno'
#base.dir = '../demo'
save.dir = './res/combine-All/'
useModel = modelLINEAR


# Genotype file name
SNP_file_name = paste0(base.dir,"/SNP.txt");
snps_location_file_name = paste0(base.dir, "/snpsloc.txt");

# Gene expression file name
expression_file_name = paste0(base.dir,'/',ct, "_GE.txt");
gene_location_file_name = paste0(base.dir, "/geneloc.txt");

# Covariates file name
# Set to character() for no covariates
covariates_file_name = paste0(base.dir,"/Covariates.txt");

# Output file name
output_file_name_cis = tempfile();
output_file_name_tra = tempfile();

pvOutputThreshold_cis = 1;
pvOutputThreshold_tra = 1e-2;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

# Distance for local gene-SNP pairs
cisDist = 1e6;

## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = " ";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = " ";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);


## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = " ";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
    cvrt$LoadFile(covariates_file_name);
}

## Run the analysis
snpspos = read.table(snps_location_file_name,
                        header = TRUE, stringsAsFactors = FALSE);
snpspos$snpid = rownames(snpspos)
snpspos = snpspos[c('snpid','chr','pos')]
genepos = read.table(gene_location_file_name,
                        header = TRUE, stringsAsFactors = FALSE);
genepos$geneid = rownames(genepos)
genepos = genepos[c('geneid','chr','left','right')]

snp_samples <- snps$columnNames
gene_samples <- gene$columnNames
common_samples <- intersect(snp_samples, gene_samples)

snps$ColumnSubsample(match(common_samples,snps$columnNames))

gene$ColumnSubsample(match(common_samples,gene$columnNames))

cvrt$ColumnSubsample(match(common_samples,cvrt$columnNames))

me = Matrix_eQTL_main(
        snps = snps,
        gene = gene,
        cvrt = cvrt,
        output_file_name  = output_file_name_tra,
        pvOutputThreshold = pvOutputThreshold_tra,
        useModel = useModel,
        errorCovariance = errorCovariance,
        verbose = TRUE,
        output_file_name.cis = output_file_name_cis,
        pvOutputThreshold.cis = pvOutputThreshold_cis,
        snpspos = snpspos,
        genepos = genepos,
        cisDist = cisDist,
        pvalue.hist = TRUE,
        min.pv.by.genesnp = TRUE,
        noFDRsaveMemory = FALSE);


saveRDS(me,paste0(save.dir,'/',ct, "_eQTL.rds"))

res = me$cis$eqtls
res = res[res$FDR < 0.05,]
write.csv(res,paste0(save.dir,'/',ct, "_cis_fdr005.csv"))
