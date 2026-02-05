import pandas as pd
import anndata as ad
import scanpy as sc
import numpy as np
# import doubletdetection
import os
import scanpy.external as sce
import scipy.stats as stats
import seaborn as sns
import statsmodels.api as sm
import warnings
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import math
import scipy
import matplotlib as mpl
import statsmodels.stats.weightstats as sw


meta = pd.read_csv('All_harmony_35_27_1.3_meta.csv', index_col=0)
raw = sc.read('All_109T1D_56HC_raw.h5ad')
raw = raw[meta.index]

# Assign cell type labels to the raw AnnData object
raw.obs['Celltype1'] = meta['Celltype1']

# Normalize raw counts (library-size normalization)
sc.pp.normalize_total(raw, target_sum=1e4)

# -------------------------
# Pseudo-bulk average by sample (all cells)
# -------------------------
res = pd.DataFrame(columns=raw.var_names, index=raw.obs['sample'].cat.categories)
for clust in raw.obs['sample'].cat.categories:
    res.loc[clust] = raw[raw.obs['sample'].isin([clust]), :].X.mean(0)
res.to_csv('./input/sc-scale/01raw/{}.csv'.format('Pan'))

# -------------------------
# Pseudo-bulk average by sample within each Celltype1
# -------------------------
for ct in raw.obs['Celltype1'].unique():
    temp = raw[raw.obs['Celltype1'] == ct]
    res = pd.DataFrame(columns=temp.var_names, index=temp.obs['sample'].cat.categories)
    for clust in temp.obs['sample'].cat.categories:
        res.loc[clust] = temp[temp.obs['sample'].isin([clust]), :].X.mean(0)
    res.to_csv('./input/sc-scale/01raw/{}.csv'.format(ct))