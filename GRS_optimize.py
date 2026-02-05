# git clone git@github.com:JassicaZ/backward_selection_GRS_tools.git
# git clone git@github.com:JassicaZ/grsopt.git
# pip install ./grsopt
# pip install ./backward_selection_GRS_tools

import grsopt
from backward_selection_GRS_tools import myutils
import pandas as pd
import numpy as np

# Load the input data
merged_df = pd.read_csv('../data/feact_sele_input_snp+dis.csv', index_col='id')
beta = pd.read_csv('../data/GRS_model_formula_input.csv')

# Extract SNP data
snp_df = merged_df[beta['loci'].tolist()]

# Multiply SNP values by their corresponding beta values
beta_dict = beta.set_index('loci')['beta'].to_dict()
for i in snp_df.columns:
    snp_df[i] = snp_df[i] * beta_dict[i]

# Perform backward elimination to select features
X = snp_df.dropna(axis=1)
y = merged_df.dropna(axis=1)['Disease']

final_selected, final_history = myutils.backward_elimination(X, y, min_features=5, cv=5)


# Load datasets
selected_matrix = merged_df[final_selected]
train_merged_df = pd.read_csv('../data/feact_sele_input_snp+dis.csv', index_col='id')
val_merged_df = pd.read_csv('../data/val_feact_sele_input_snp+dis.csv', index_col="id")
beta_orign = pd.read_csv('../data/GRS_model_formula_input.csv', index_col='loci')  # Original beta values

# Extract feature matrices for training and validation
X_val = val_merged_df.iloc[:, 2:]
X_train = selected_matrix

# Retain only the intersecting loci between training and validation sets
geno_interect = list(set(X_val.columns) & set(X_train.columns))
X_val = X_val[geno_interect]

# Ensure the loci are in the same order in both datasets
X_val = X_val[X_train.columns]

# Replace NaN values in the validation set with 0 (no contribution to calculations)
X_val = X_val.fillna(0)

# Map disease labels to numerical values
mapping = {"HC": 0, "T1D": 1}
y_train = train_merged_df['Disease'].map(mapping)
y_val = val_merged_df['Disease'].map(mapping)

# Align the order of beta values with the training feature matrix
beta_orign = beta_orign.loc[X_train.columns]

# Optimize beta values
grsopter = grsopt.grsoptimizer(original_betas=beta_orign['beta'],
                               regularization_strength=0.1)
result = grsopter.optimize(X_train.reset_index(drop=True), y_train)

# Calculate GRS for the validation set
GRS_val_ori = np.dot(X_val, beta_orign['beta'])
GRS_val_opt = np.dot(X_val, result['optimized_betas'])

