import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from grsopt import calculate_auc
from sklearn.metrics import roc_auc_score, average_precision_score, roc_curve


grs_val_all = pd.read_csv("./grs_input/grs_val_all-LD.csv")
grs_train_all = pd.read_csv("./grs_input/grs_train_all-LD.csv")

# ROC for validation cohort
mpl.rcParams["font.family"] = "Arial"
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42

tmp = grs_val_all[["GRS-us", "Disease"]].dropna()

y_true = tmp["Disease"].map({"HC": 0, "T1D": 1})
y_score = tmp["GRS-us"]

fpr, tpr, _ = roc_curve(y_true, y_score)
auc = roc_auc_score(y_true, y_score)

plt.figure(figsize=(6, 5))

plt.plot(
    fpr,
    tpr,
    label=f"GRS-us (AUC = {auc:.2f})"
)

plt.plot(
    [0, 1],
    [0, 1],
    linestyle="--"
)

plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("ROC Curve of GRS-us")
plt.legend()
plt.tight_layout()

#plt.show()
plt.savefig(
    "./grs_output/GRS_performance_val-ld.pdf",
    format="pdf",
    bbox_inches="tight"
)

#ROC for discovery cohort

tmp = grs_train_all[["GRS-us", "Disease"]].dropna()

y_true = tmp["Disease"].map({"HC": 0, "T1D": 1})
y_score = tmp["GRS-us"]

fpr, tpr, _ = roc_curve(y_true, y_score)
auc = roc_auc_score(y_true, y_score)

plt.figure(figsize=(6, 5))

plt.plot(
    fpr,
    tpr,
    label=f"GRS-us (AUC = {auc:.2f})"
)

plt.plot(
    [0, 1],
    [0, 1],
    linestyle="--"
)

plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("ROC Curve of GRS-us")
plt.legend()
plt.tight_layout()

#plt.show()
plt.savefig(
    "./grs_output/GRS_performance_train-LD.pdf",
    format="pdf",
    bbox_inches="tight"
)


## Feature importance

def snp_importance(X, y, betas):
    cols = list(X.columns)
    full = calculate_auc(X, y, betas, cols)      # AUC using all SNPs
    rows = []
    for snp in cols:
        sub = [c for c in cols if c != snp]      # Remove the current SNP
        auc_wo = calculate_auc(X, y, betas, sub)
        rows.append({"SNP": snp, "score_with": full,
                     "score_without": auc_wo,
                     "importance": full - auc_wo})
    return (pd.DataFrame(rows)
              .set_index("SNP")
              .sort_values("importance", ascending=False))

train_df=pd.read_csv('./grs_output/Discovery_risk_matrix-selected.csv', index_col=0)
df= pd.read_csv('./resgrs_outputult/optimized_beta-LD.csv')
X_train= train_df.T

train_cli=pd.read_csv('./grs_output/Dis_cli.csv', index_col=0)
y_train=train_cli.loc[train_df.columns]['Disease']
mapping = {"HC": 0, "T1D": 1}
y_train =y_train.map(mapping)

imp_train = snp_importance(X_train, y_train, df['beta_opt'])

import pandas as pd
import matplotlib.pyplot as plt

mpl.rcParams["font.family"] = "Arial"
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42


df_sorted = imp_train.sort_values("importance", ascending=False).copy()

# Compute relative importance
max_importance = df_sorted["importance"].iloc[0]
df_sorted["relative_importance"] = (
    df_sorted["importance"] / max_importance
)

# Top 20 SNPs
top20 = df_sorted.head(20)

# Horizontal bar plot
plt.figure(figsize=(6, 7))

plt.barh(
    top20.index,
    top20["relative_importance"],
    #height=0.6
)

# Place the most important SNP at the top
plt.gca().invert_yaxis()

plt.xlabel("Relative importance")
plt.ylabel("SNP")
plt.title("Top 10 SNPs by relative importance")

plt.tight_layout()
plt.savefig(
    "./grs_output/importance_train-LD-20.pdf",
    format="pdf",
    bbox_inches="tight"
)
plt.show()





