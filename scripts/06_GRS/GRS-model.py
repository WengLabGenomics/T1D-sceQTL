from grsopt import backward_selection, SelectionConfig, OptimizerConfig, forward_selection
import pandas as pd
import numpy as np

# 读取数据
train_df=pd.read_csv('./grs_input/us_only/Discovery_risk_matrix-us-LD.csv', index_col=0)
val_df=pd.read_csv('./grs_input/us_only/Validation_risk_matrix-us-LD.csv',index_col=0)

beta_df =pd.read_csv('./grs_input/us_only/SNP_beta-us-LD.csv',index_col='SNP')

snp_orgin= list(set(train_df.index) & set(val_df.index) & set(beta_df.index))
X_train= train_df.loc[snp_orgin,].T
X_val= val_df.loc[snp_orgin,].T
beta_orign=beta_df.loc[snp_orgin,]

train_cli=pd.read_csv('./grs_input/Dis_cli.csv', index_col=0)
val_cli=pd.read_csv('./grs_input/Val_cli.csv', index_col=0)

y_train=train_cli.loc[train_df.columns]['Disease']
y_val=val_cli.loc[val_df.columns]['Disease']

mapping = {"HC": 0, "T1D": 1}
y_train =y_train.map(mapping)
y_val =y_val.map(mapping)

# Feature selection
from grsopt import grsoptimizer

sel = SelectionConfig(
    direction="backward", 
    cv_folds=5,
    start_n=None,           
    end_n=10,               
    margin=0.0,
    retention_threshold=0.8, 
)
opt = OptimizerConfig(max_iter=200, n_initial=3, regularization_strength=0.05,
                      noise_scale=0.8, loss='rank')

res = backward_selection(X_train, y_train, beta_orign['b'], sel, opt)

final_selected = res["selected_snps"]
X_train_selected= X_train[final_selected]
beta_orign_selected = beta_orign.loc[final_selected]

final_selected = res["selected_snps"]
X_train_selected= X_train[final_selected]
beta_orign_selected = beta_orign.loc[final_selected]

# Formal beta optimization
opt = grsoptimizer(
    original_betas=beta_orign_selected['b'],
    regularization_strength=0.1,
    max_iter=200, n_initial=5, noise_scale=0.01,
    loss='rank',
    verbose=True
)
result = opt.optimize(X_train_selected, y_train)

df = pd.DataFrame(
    result['optimized_betas'],
      index=final_selected,
   columns=['beta_opt']
)

grs_train = np.dot(X_train_selected.to_numpy(), df['beta_opt'].to_numpy())
grs_val = np.dot(X_val[final_selected].to_numpy(), df['beta_opt'].to_numpy())

grs_train_df = pd.DataFrame({'GRS-us': grs_train, 'Disease': y_train})
grs_val_df = pd.DataFrame({'GRS-us': grs_val, 'Disease': y_val})

df.to_csv('./grs_output/optimized_beta-LD.csv')