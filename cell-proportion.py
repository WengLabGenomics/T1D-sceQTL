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


adata = sc.read('All_harmony_35_27_1.3_withcli.h5ad.h5ad')

def formatter(x, pos):
    """The two args are the value and tick position."""
    if ((x > 0.9) or (x == 0)):
        return '%1.0f' % x
    else:
        res = '%1.1f' % x
        return res


def convert_pvalue_to_asterisks(pvalue):
    if pvalue <= 0.0001:
        return "****"
    elif pvalue <= 0.001:
        return "***"
    elif pvalue <= 0.01:
        return "**"
    elif pvalue <= 0.05:
        return "*"
    return "ns"


def get_proportion_stats(ind_perc, group1, group2, covar, ct_name):
    """
    Weighted least squares:
      ct_perc ~ 1 + I(group2)
    weights = ind_count_sum
    Multiple testing correction: Bonferroni across cell types.
    """
    print('{0} versus {1}'.format(group1, group2))
    all_out = pd.DataFrame()
    ind_perc0 = ind_perc[ind_perc[covar].isin(group1 + group2)]

    Effect = {}
    Pval = {}

    for ct_i in list(range(len(ind_perc0[ct_name].cat.categories))):
        ct = ind_perc0[ct_name].cat.categories[ct_i]

        ct_perc = ind_perc0.ct_perc[ind_perc0[ct_name] == ct]
        weights = ind_perc0.ind_count_sum[ind_perc0[ct_name] == ct]

        disease = ind_perc0[covar][ind_perc0[ct_name] == ct].astype("str")
        disease.values[disease.isin(group1)] = 0
        disease.values[disease.isin(group2)] = 1

        disease = sm.add_constant(disease.astype(float))
        est = sm.WLS(ct_perc.astype(float), disease, weights=weights)
        est = est.fit()

        Effect[ct] = est.params[1]
        Pval[ct] = est.pvalues[1]

        all_out = pd.concat(
            [all_out,
             pd.DataFrame({"Cell": str(ct),
                           "Beta": str(Effect[ct]),
                           "Pval": str(Pval[ct]),
                           "sig": convert_pvalue_to_asterisks(Pval[ct])},
                          index=[ct])],
            ignore_index=True
        )

    all_out['Pval'] = all_out['Pval'].astype(float)

    # Bonferroni correction across tested cell types
    m = len(all_out['Pval'])
    all_out['Pval_bonf'] = np.minimum(all_out['Pval'] * m, 1.0)
    all_out['sig_bonf'] = all_out['Pval_bonf'].apply(convert_pvalue_to_asterisks)

    # Return raw p-value dict (as you used later) + also bonferroni dict for convenience
    p_out = all_out.set_index('Cell')['Pval'].to_dict()
    p_bonf_out = all_out.set_index('Cell')['Pval_bonf'].to_dict()

    # If you still want to "highlight" significant ones (Bonferroni)
    def highlight_sig_bonf(p):
        sig = p < (0.05 / len(p))
        return ['background-color: yellow' if v else '' for v in sig]

    try:
        display(all_out.style.apply(highlight_sig_bonf, subset=['Pval']))
    except NameError:
        pass

    return Effect, p_out, p_bonf_out, all_out


# ---- orders / colors (keep your original definitions) ----
mye_order = ['CD14Mono','CD16Mono','cDC1','cDC2','pDC']
lym_order = ['CD4T_Naive','CD4T_Mem','CD4T_reg',
             'CD8T_Naive','CD8T_Mem','CD8T_MAIT','CD8T_Cytotoxic',
             'NK_Dim','NK_Bright']
B_order = ['B_Naive','B_Mem','B_Atypical','B_Plasma','B_Plasma_cycling']
M_order  = mye_order + lym_order + B_order + ['T_Prolif']

mye_color = ['#FFDAB9','#CD6600', '#8B4500', '#E9967A','#FF8C00']
lym_color = ['#0000EE','#00008B', '#1E90FF',
             '#00F5FF' ,'#98F5FF', '#00CED1', '#7AC5CD',
             '#AB82FF','#E066FF']
B_color = ['#8B0000',  '#FF69B4', '#CD5C5C',  '#FF0000',  '#FFB6C1']
M_color = mye_color + lym_color + B_color + ['#696969']

MasterORDER_dict = {
    'Celltype2': M_order
}
colorrSS_dict = {
    'Celltype2': M_color
}


def Celltype_pro_plot(covar, ct_name, order_list, group1_number, group2_number):
    MasterORDER = MasterORDER_dict[ct_name]
    colorrSS = colorrSS_dict[ct_name]
    col_wrap_dict = {'Celltype2': 5}
    colorrs_desat = sns.color_palette(colorrSS, desat=1)

    # Use adata (fix: adata2 was undefined)
    adata_obs_small = adata.obs.copy()
    adata_obs_small[covar] = adata_obs_small[covar].astype('str').astype('category')

    ind_count = adata_obs_small.groupby(['sample', covar, ct_name], observed=True)[ct_name].count()
    ind_count_sums = ind_count.groupby(level=[0]).sum().reset_index(name="counts")

    ind_perc = (ind_count / ind_count.groupby(level=[0]).transform(sum)) * 100
    ind_perc = ind_perc.reset_index(name="ct_perc")

    # Add weights to WLS
    ind_perc['counts'] = ind_count.values.tolist()
    ind_perc['ind_count_sum'] = list(np.zeros(len(ind_count.values.tolist()), dtype=int))

    # Add total sums per individual
    for ii in range(len(ind_count_sums)):
        ind_perc.loc[ind_perc['sample'] == ind_count_sums['sample'][ii], 'ind_count_sum'] = ind_count_sums.counts[ii]

    ind_perc[ct_name] = ind_perc[ct_name].astype('category')
    ind_perc[covar] = ind_perc[covar].astype('str').astype('category')

    perc_plot = sns.catplot(
        x=covar, y='ct_perc',
        hue=ct_name, data=ind_perc, kind='violin',
        col_order=MasterORDER, order=order_list,
        col=ct_name, col_wrap=col_wrap_dict[ct_name],
        cut=0, dodge=False,
        height=2, aspect=0.8, linewidth=0.5,
        sharex=False, sharey=False, palette=colorrSS
    )

    perc_plot.set_yticklabels(size=10)
    perc_plot.set_xticklabels(size=10)

    for ct_i in list(range(len(MasterORDER))):
        ct = MasterORDER[ct_i]
        sns.stripplot(
            x=covar, y="ct_perc",
            data=ind_perc[ind_perc[ct_name] == ct],
            order=order_list,
            color=colorrs_desat[ct_i],
            jitter=True, size=2, linewidth=0.5, ax=perc_plot.axes[ct_i]
        )
        perc_plot.axes[ct_i].get_yaxis().label.set_visible(False)
        perc_plot.axes[ct_i].get_xaxis().label.set_visible(False)
        perc_plot.axes[ct_i].tick_params(axis="y", direction='in', length=3, width=1, pad=2)
        perc_plot.axes[ct_i].tick_params(axis="x", direction='out', length=3, width=1, pad=1)
        perc_plot.axes[ct_i].yaxis.set_major_formatter(mtick.FuncFormatter(formatter))
        perc_plot.axes[ct_i].set_title(ct)
        perc_plot.fig.tight_layout()

    perc_plot.set_xticklabels(rotation=90)
    perc_plot.fig.subplots_adjust(wspace=0.25, hspace=0.9)
    perc_plot.set_xticklabels(order_list)

    # WLS stats + Bonferroni
    HLA_effect, HLA_Pval, HLA_Pval_bonf, out_df = get_proportion_stats(
        ind_perc,
        group1=[order_list[group1_number]],
        group2=[order_list[group2_number]],
        covar=covar,
        ct_name=ct_name
    )
    return out_df

# Run
Celltype_pro_plot('Disease', 'Celltype2', ['HC', 'T1D'], 0, 1)