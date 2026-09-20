import os
import warnings

import numpy as np
import pandas as pd
import scanpy as sc

from sklearn.linear_model import LinearRegression


adata = sc.read('All_harmony_35_27_1.3_withcli.h5ad')

adata2 = adata[
    adata.obs['Disease'] == 'T1D'
].copy()


adata2.obs['Celltype1'] = adata2.obs['Celltype1'].astype(str)
adata2.obs['Celltype2'] = adata2.obs['Celltype2'].astype(str)


ct = 'Celltype2'   # change to 'Celltype1'

min_n_cells_total = 50   # minimum total cells for each cell type
min_n_cells = 10         # minimum cells per donor
min_n_samples = 2       

n_pcs = 20


covariates = [
    'onset_age',
    'age',
    'Duration',
    'sex',
    'HbA1c',
    'ZnT8_level',
    'IA2_level',
    'GAD_level'
]


# Continuous covariates
contin = [
    'onset_age',
    'age',
    'Duration',
    'HbA1c',
    'ZnT8_level',
    'IA2_level',
    'GAD_level'
]


dir_results = "./results/variance_explained_by_covariates"

# Create output directories
os.makedirs(
    os.path.join(
        dir_results,
        "variance_explained_fractions"
    ),
    exist_ok=True
)

os.makedirs(
    os.path.join(
        dir_results,
        "samples_included"
    ),
    exist_ok=True
)

overwrite = True

int_types = ["T1D_ct2"]

samples_included = {}

cts_to_skip = []


for int_type in int_types:

    samples_included[int_type] = {}

    # cell types + whole atlas
    subsets = (
        sorted(
            adata2.obs[ct]
            .dropna()
            .unique()
        )
        + ["whole_atlas"]
    )


    for subset in subsets:

        subset_no_space = subset.replace(" ", "_")


        print(
            f"\nWorking on {int_type}, {subset}..."
        )

        output_file = os.path.join(
            dir_results,
            "variance_explained_fractions",
            f"variance_explained_fractions_{subset_no_space}_{int_type}.csv"
        )

        if (not overwrite) and os.path.isfile(output_file):

            print(
                f"{output_file} already exists. Skipping."
            )

            continue

        # Initialize sample tracking


        samples_included[int_type][subset] = pd.DataFrame(
            False,
            index=adata2.obs["sample"].unique(),
            columns=covariates,
            dtype=bool
        )
        
        # Select cells

        if subset == "whole_atlas":
            subadata = adata2.copy()

            verbose = True

        else:

            subadata = adata2[
                adata2.obs[ct] == subset,
                :
            ].copy()

            verbose = False


        # minimum total cells
        if subadata.n_obs < min_n_cells_total:

            print(
                f"{subset} has fewer than "
                f"{min_n_cells_total} cells. Skipping."
            )

            cts_to_skip.append(subset)

            continue

        # PCA

        # Make sure HVG information exists
        if 'highly_variable' not in subadata.var.columns:

            raise ValueError(
                "'highly_variable' is not present in adata.var. "
                "Please verify that HVGs were defined before this analysis."
            )


        n_hvg = int(
            subadata.var['highly_variable'].sum()
        )


        if n_hvg < 2:

            print(
                f"{subset} has fewer than 2 highly variable genes. "
                "Skipping."
            )

            cts_to_skip.append(subset)

            continue

        n_comps_use = min(
            n_pcs,
            subadata.n_obs - 1,
            n_hvg - 1
        )


        if n_comps_use < 1:

            print(
                f"Not enough dimensions for PCA in {subset}. "
                "Skipping."
            )

            continue


        sc.tl.pca(
            subadata,
            n_comps=n_comps_use,
            use_highly_variable=True
        )


        emb_name = "X_pca"

        n_comps = subadata.obsm[
            emb_name
        ].shape[1]

        # Create cell-level table

        comp_sample_df = pd.DataFrame(
            index=subadata.obs_names
        )


        comp_sample_df["sample"] = (
            subadata.obs["sample"].values
        )


        comp_sample_df["n_cells"] = 1


        # PC scores
        for comp in range(n_comps):

            comp_sample_df[
                f"comp{comp}"
            ] = subadata.obsm[
                emb_name
            ][:, comp]


        # Clinical variables
        for cov in covariates:

            if cov not in subadata.obs.columns:

                raise KeyError(
                    f"{cov} is not present in adata.obs"
                )

            comp_sample_df[cov] = (
                subadata.obs[cov].values
            )


        for cov in covariates:

            n_unique = (
                comp_sample_df
                .groupby(
                    "sample",
                    observed=True
                )[cov]
                .nunique(
                    dropna=True
                )
            )


            problematic = n_unique[
                n_unique > 1
            ]


            if len(problematic) > 0:

                warnings.warn(
                    f"{cov}: {len(problematic)} samples have "
                    "more than one value within the same donor."
                )

        # Aggregate to donor level

        agg_dict = {
            "n_cells": "sum"
        }


        for comp in range(n_comps):

            agg_dict[
                f"comp{comp}"
            ] = "mean"


        for cov in covariates:

            agg_dict[cov] = "first"


        sample_df = (
            comp_sample_df
            .groupby(
                "sample",
                observed=True
            )
            .agg(agg_dict)
        )


        # donor must contain at least min_n_cells
        sample_df = sample_df.loc[
            sample_df["n_cells"] >= min_n_cells
        ].copy()


        print(
            f"{subset}: "
            f"{sample_df.shape[0]} donors after "
            f"minimum-cell filtering."
        )


        if sample_df.shape[0] < min_n_samples:

            print(
                f"Only {sample_df.shape[0]} samples available "
                f"for {subset}. Skipping."
            )

            cts_to_skip.append(subset)

            continue

        # Variance explained

        var_explained = pd.DataFrame(
            np.nan,
            index=range(n_comps),
            columns=covariates,
            dtype=float
        )


        total_variance = pd.DataFrame(
            np.nan,
            index=range(n_comps),
            columns=covariates,
            dtype=float
        )


        # Loop over covariates

        for cov in covariates:

            # Continuous covariate

            if cov in contin:

                x_series = pd.to_numeric(
                    sample_df[cov],
                    errors="coerce"
                )


                valid = x_series.notna()


                if valid.sum() < min_n_samples:

                    print(
                        f"{subset} - {cov}: only "
                        f"{valid.sum()} valid samples. Skipping."
                    )

                    continue


                x = (
                    x_series
                    .loc[valid]
                    .to_numpy(dtype=float)
                    .reshape(-1, 1)
                )


                if verbose:

                    print(
                        f"Treating {cov} as continuous "
                        f"(n={valid.sum()})"
                    )


            # Categorical covariate

            else:

                x_series = (
                    sample_df[cov]
                    .astype("string")
                    .replace(
                        {
                            "nan": pd.NA,
                            "NaN": pd.NA,
                            "None": pd.NA,
                            "none": pd.NA,
                            "NA": pd.NA,
                            "N/A": pd.NA,
                            "": pd.NA
                        }
                    )
                )


                valid = x_series.notna()


                if valid.sum() < min_n_samples:

                    print(
                        f"{subset} - {cov}: only "
                        f"{valid.sum()} valid samples. Skipping."
                    )

                    continue


                x_valid = x_series.loc[
                    valid
                ]


                # Must contain >=2 categories
                if x_valid.nunique() < 2:

                    print(
                        f"{subset} - {cov}: only one category "
                        "is present. Skipping."
                    )

                    continue


                x = pd.get_dummies(
                    x_valid,
                    drop_first=True,
                    dtype=float
                ).to_numpy()


                if verbose:

                    print(
                        f"Treating {cov} as categorical "
                        f"(n={valid.sum()}, "
                        f"categories={x_valid.nunique()})"
                    )


            # Store included donors

            valid_samples = sample_df.index[
                valid
            ]


            samples_included[
                int_type
            ][subset][cov] = (
                samples_included[
                    int_type
                ][subset]
                .index
                .isin(valid_samples)
            )

            # Loop over PCs

            for comp in range(n_comps):

                y_true = (
                    sample_df
                    .loc[
                        valid,
                        f"comp{comp}"
                    ]
                    .to_numpy(
                        dtype=float
                    )
                )


                model = LinearRegression(
                    fit_intercept=True
                )


                model.fit(
                    x,
                    y_true
                )


                y_pred = model.predict(
                    x
                )


                # variance predicted by this covariate
                var_explained.loc[
                    comp,
                    cov
                ] = np.var(
                    y_pred
                )


                total_variance.loc[
                    comp,
                    cov
                ] = np.var(
                    y_true
                )

        # Sum variance across PCs

        total_variance_explained = (
            var_explained.sum(
                axis=0,
                skipna=True,
                min_count=1
            )
        )


        total_variance_observed = (
            total_variance.sum(
                axis=0,
                skipna=True,
                min_count=1
            )
        )


        # avoid division by zero
        total_variance_observed = (
            total_variance_observed
            .replace(
                0,
                np.nan
            )
        )

        total_variance_explained_fractions = (
            total_variance_explained
            /
            total_variance_observed
        )


        total_variance_explained_fractions = (
            total_variance_explained_fractions
            .sort_values(
                ascending=False
            )
        )


        total_variance_explained_fractions.name = (
            "variance_explained_fraction"
        )


        # Save

        total_variance_explained_fractions.to_csv(
            output_file
        )


        samples_included[
            int_type
        ][subset].to_csv(
            os.path.join(
                dir_results,
                "samples_included",
                f"samples_included_{subset_no_space}_{int_type}.csv"
            )
        )


        print(
            f"Finished: {subset}"
        )


# Summary

print("\nAnalysis completed.")

if len(cts_to_skip) > 0:

    print(
        "Skipped subsets:",
        sorted(set(cts_to_skip))
    )