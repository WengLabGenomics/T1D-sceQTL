import pandas as pd
import anndata as ad
import scanpy as sc
import numpy as np
import doubletdetection
import os
import scanpy.external as sce

maindir = os.path.expanduser('~/workspace/T1D/cellranger_out/')
date = 'all'
subdir = 'scanpy_process_data'

match = pd.read_csv(maindir + '/scanpy_useid.csv')
mix_id = match['mix_number'].to_list()

clf = doubletdetection.BoostClassifier(
    n_iters=10,
    clustering_algorithm="louvain",
    standard_scaling=True,
    pseudocount=0.1,
    n_jobs=-1,
)

for i in mix_id:
    company_id = list(match.loc[match['mix_number'] == i, 'company_number'])

    for j in company_id:
        # Read 10x matrix for each sample
        temp = sc.read_10x_mtx(
            '{}/{}/outs/per_sample_outs/{}/count/sample_filtered_feature_bc_matrix/'.format(maindir, i, j)
        )
        temp.obs['sample'] = str(j)
        temp.obs['batch'] = str(i)
        temp.obs_names = list(str(j) + '_' + temp.obs_names)

        # Doublet detection
        doublets = clf.fit(temp.X).predict(p_thresh=1e-16, voter_thresh=0.5)
        doublet_score = clf.doublet_score()
        temp.obs["doublet"] = doublets
        temp.obs["doublet_score"] = doublet_score

        if j == company_id[0]:
            adata = temp
        else:
            adata = ad.concat([adata, temp], join='outer')

# Remove doublets
adata = adata[adata.obs['doublet'] == 0].copy()

file_path = maindir + '/' + subdir
if os.path.exists(file_path):
    print('exist')
else:
    os.mkdir(file_path)

adata.write(file_path + '/T1D_merge_outer.h5ad')

# QC metrics
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate mitochondrial genes
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

sc.pl.violin(
    adata,
    ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
    jitter=0,
    multi_panel=True,
    groupby='sample',
    save='QC.pdf'
)

adata = adata[adata.obs.n_genes_by_counts > 200, :].copy()
# adata = adata[adata.obs.n_genes_by_counts < 7000, :].copy()
adata = adata[adata.obs.pct_counts_mt < 15, :].copy()

# Remove mitochondrial/ribosomal/heat-shock/IEG-related genes (based on provided list)
IEG = pd.read_table(maindir + '/scanpy_process_data/IEG.txt', header=None)
adata = adata[:, adata.var_names.isin(list(set(adata.var_names) - set(IEG[0])))].copy()

# Normalization & HVG
sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.log1p(adata)
adata.raw = adata
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')

# Harmony integration
# NOTE: key should be a string column name, not a list
sce.pp.harmony_integrate(adata, key='batch')
print('harmony is done!')

adata.write(file_path + '/{}_harmony_all.h5ad'.format(date))

n_n = 35
n_pc = 30
sc.pp.neighbors(adata, n_neighbors=n_n, n_pcs=n_pc, use_rep='X_pca_harmony')
sc.tl.umap(adata)

sc.settings.verbosity = 3  # errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

sc.pl.umap(adata, color=['sample', 'batch'], save='inte_batch.pdf')

res = 1.3
sc.tl.leiden(adata, resolution=res)
sc.pl.umap(adata, color=['leiden'], save='neig{}_pc{}_leiden_{}.pdf'.format(n_n, n_pc, res))

sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon', use_raw=True)
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, save='res{}_leiden_deg.pdf'.format(res))

adata.write(file_path + '/{}_harmony_{}_{}_{}_withcli.h5ad'.format(date, n_n, n_pc, res))