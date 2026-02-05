import pandas as pd
import anndata as ad
import scanpy as sc
import numpy as np
import doubletdetection
import os
import scanpy.external as sce
import scipy.stats as stats
import statsmodels.api as sm
import warnings
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import math
import scipy
import statsmodels.stats.weightstats as sw
#import utils
from sklearn.linear_model import LinearRegression
from scipy.stats import pearsonr
from scipy.stats import entropy

# --- load data ---
adata = sc.read('All_harmony_35_27_1.3_withcli.h5ad')  # fix: remove duplicated .h5ad
adata.obs = adata.obs.replace('99999', pd.NA, regex=True)

# --- basic recoding ---
adata.obs['GAD_state'] = adata.obs['GAD_state'].astype('str')
adata.obs['IA2_state'] = adata.obs['IA2_state'].astype('str')
adata.obs['sex'] = adata.obs['sex'].astype('str')
adata.obs['HLA_type_new'] = adata.obs['HLA_type_new'].astype('str')

adata.obs.loc[adata.obs['GAD_state'] == 'Low', 'GAD_state'] = 0
adata.obs.loc[adata.obs['GAD_state'] == 'High', 'GAD_state'] = 1
adata.obs.loc[adata.obs['IA2_state'] == 'Low', 'IA2_state'] = 0
adata.obs.loc[adata.obs['IA2_state'] == 'High', 'IA2_state'] = 1
adata.obs.loc[adata.obs['HLA_type_new'] == 'S', 'HLA_type_new'] = 1
adata.obs.loc[adata.obs['HLA_type_new'] == 'N', 'HLA_type_new'] = 0

# --- Ab count ---
adata.obs[['IA2_ID','GAD_ID','ZnT8_ID']] = adata.obs[['IA2_ID','GAD_ID','ZnT8_ID']].astype(str)
tmp = adata.obs[['IA2_ID','GAD_ID','ZnT8_ID']].copy()
tmp.replace('nan', np.nan, inplace=True)
tmp.fillna(0, inplace=True)
tmp = tmp.astype(int)
adata.obs['Ab'] = tmp.sum(axis=1).astype(int)

adata.obs['Disease'] = adata.obs['Disease'].astype('str')
adata2 = adata[adata.obs['Disease'] == 'T1D'].copy()
adata2.obs['Celltype2'] = adata2.obs['Celltype2'].astype(str)
adata2.obs['Celltype1'] = adata2.obs['Celltype1'].astype(str)

ct = 'Celltype2'  # or 'Celltype1'

min_n_cells_total = 50  # total cells per subset
min_n_cells = 10        # cells per sample
min_n_samples = 2

covariates = [
    'onset_age','age','Duration','sex','BMI','HbA1c','ZnT8_level','IA2_level','GAD_level',
    'ZnT8_ID','IA2_ID','GAD_ID','WHR','Ab','IA2_state','GAD_state','HLA_type_new'
]
contin = ['onset_age','age','Duration','BMI','HbA1c','ZnT8_level','IA2_level','GAD_level','WHR','GRS_Disease','GRS_Age']

dir_results = "./results/variance_explained_by_covariates"
n_pcs = 20
int_types = ["T1D_ct2"]

# Ensure output folders exist
os.makedirs(os.path.join(dir_results, "variance_explained_fractions"), exist_ok=True)
os.makedirs(os.path.join(dir_results, "samples_included"), exist_ok=True)

samples_included = dict()

for int_type in int_types:
    samples_included[int_type] = dict()
    # loop through all annotations (= cell types), and also include the entire
    # atlas as a "cell type":
    for subset in sorted(adata2.obs[ct].unique()) + ["whole_atlas"]:
        subset_no_space = subset.replace(" ", "_")
        samples_included[int_type][subset] = pd.DataFrame(
            index=adata.obs["sample"].unique(), columns=covariates
        )
        if not os.path.isfile(
            os.path.join(
                dir_results,
                f"variance_explained_fractions/variance_explained_fractions_{subset_no_space}_{int_type}.csv",
            )
        ):
            print(f"Working on {int_type}, {subset}...")
            # select the correct cells:
            if subset == "whole_atlas":
                subadata = adata.copy()
                verbose = True
            elif subset not in adata2.obs[ct].unique():
                raise ValueError(
                    "subset should be set either to 'Whole atlas' or to a category in your manual_ann grouped obs variable!"
                )
            else:
                subadata = adata2[adata2.obs[ct] == subset, :].copy()
                verbose = False
            if subadata.n_obs < min_n_cells_total:
                print(f"{subset} has fewer than {min_n_cells_total} cells! Skipping.")
                continue
            # select the right embedding:
            emb_name = "X_pca"
            sc.tl.pca(subadata, n_comps=n_pcs, use_highly_variable=True)
            n_comps = subadata.obsm[emb_name].shape[1]
            var_explained = pd.DataFrame(
                index=range(n_comps), columns=covariates + ["overall"]
            )
            # initiate a dataframe in which we will store the data
            # for our linear regression (i.e. the PC/latent components, + covariates).
            # Rows are cells, but we will collapse this to samples below
            comp_sample_df = pd.DataFrame(index=subadata.obs.index)
            comp_sample_df["sample"] = subadata.obs["sample"]
            agg_dict = {"sample": "count"}  # this will be number of cells
            for comp in range(n_comps):
                # store component scores per cell
                comp_sample_df[f"comp{comp}"] = subadata.obsm[emb_name][:, comp]
                # we will aggregate these later by taking the mean per sample
                agg_dict[f"comp{comp}"] = "mean"
            for cov in covariates:
                if cov in ["log10_total_counts", "mito_frac"]:
                    # store values
                    comp_sample_df[cov] = subadata.obs[cov]
                    # we will aggregate by taking the mean
                    agg_dict[cov] = "mean"
                else:
                    # for all other covariates: these are sample-level
                    # covariates, so we will take the "first" overservation
                    # in the sample (which should be the only)
                    comp_sample_df[cov] = subadata.obs[cov]
                    agg_dict[cov] = "first"
            sample_df = (
                comp_sample_df.groupby("sample")
                .agg(agg_dict)
                .rename(columns={"sample": "n_cells"})
            )
            sample_df = sample_df.loc[
                sample_df.n_cells >= min_n_cells,
            ].copy()
            if sample_df.shape[0] < min_n_samples:
                print(
                    f"Only {sample_df.shape[0]} samples available for {subset}. Skipping."
                )
                cts_to_skip.append(subset)
                continue
            for comp in range(n_comps):
                # store the component values (for all samples i.e. unfiltered)
                y_true_unfiltered = sample_df.loc[:, f"comp{comp}"].values
                # and store variance of y_true as "overall" variance
                var_explained.loc[f"comp{comp}", "overall"] = np.var(y_true_unfiltered)
                # and the covariate as fixed variable
                for cov in covariates:
                    # store covariate observations under x
                    sample_df[cov] = sample_df[cov].apply(pd.to_numeric, errors='coerce')
                    x = sample_df[cov].values.copy()
                    # store samples to which they match
                    x_samples = sample_df.index
                    # check which of these samples have no observation (e.g.
                    # because BMI was unknown, or age, etc.)
                    # (the function used below checks for different kinds of
                    # nas, e.g. np.nan, "nan", None, "None" etc.)
                    x = x.astype(float)
                    x_nans = np.isnan(x)
                    # now keep only xs that have real observations
                    x = x[~x_nans]
                    if len(x) < 2:
                        continue
                    # filter samples according to x filtering
                    x_samples = x_samples[~x_nans]
                    # and store which samples were included in our samples_included
                    # dictionary, for later reference (this is our "n")
                    samples_included[int_type][subset][cov] = samples_included[
                        int_type
                    ][subset].index.isin(x_samples.tolist())
                    # filter y_true according to x's filtering
                    y_true = y_true_unfiltered[~x_nans].reshape(-1, 1)
                    #if x.dtype in ["float32", "float", "float64"]:
                    if cov in contin:
                        x = x.reshape(-1, 1)
                        # print that we are treating as numerical (only for first comp,
                        # so that we don't print the same thing many times)
                        if comp == 0 and verbose:
                            print(f"treating {cov} as continuous variable")
                    # otherwise we are dealing with a categorical...
                    else:
                        # if it has only one category, there is 0 variance and
                        # we cannot perform linear regression. In that case,
                        # move on to the next covariate.
                        if len(set(x)) == 1:
                            var_explained.loc[comp, cov] = np.nan
                            continue
                        # Otherwise, convert x to dummied variable:
                        # print that we are converting to dummy
                        # (only do it for the first comp, otherwise we print the same thing
                        # many times)
                        if comp == 0 and verbose:
                            print(f"converting {cov} to dummy variable")
                        # drop_first means we ensure that are x is full rank,
                        # and we only encode all-1 categories
                        x = pd.get_dummies(x, drop_first=True)
                    lrf = LinearRegression(fit_intercept=True).fit(
                        x,
                        y_true,
                    )
                    # predict y based on the fit linear model
                    y_pred = lrf.predict(x)
                    # and store the variance of the predicted y, this is the
                    # "variance explained" by the covariate, for this component
                    var_explained.loc[comp, cov] = np.var(y_pred)
                total_variance_explained = np.sum(var_explained, axis=0).sort_values(
                ascending=False
            )
            # divide this by the total variance that was observed in the
            # components, to get fraction of variance explained
            total_variance_explained_fractions = (
                total_variance_explained / total_variance_explained["overall"]
            )
            total_variance_explained_fractions.to_csv(
                os.path.join(
                    dir_results,
                    f"variance_explained_fractions/variance_explained_fractions_{subset_no_space}_{int_type}.csv",
                )
            )
            samples_included[int_type][subset].to_csv(
                os.path.join(
                    dir_results,
                    f"samples_included/samples_included_{subset_no_space}.csv",
                )
            )