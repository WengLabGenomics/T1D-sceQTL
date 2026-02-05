library(dplyr)
library(data.table)
library(stringr)

## First, check how many genes remain after QC for each cell type
for (ct in c('Pan','CD8T', 'NK', 'CD4T', 'CD14Mono', 'B', 'CD16Mono', 'pDC', 'cDC','T_Prolif', 'Plasmablast')){
  data = data.frame(fread(paste0('./input/sc-scale/01raw/', ct, '.csv'), header = FALSE))
  rownames(data) <- data$V1
  data = data[, 2:ncol(data)]
  keep = colSums(data > 0) >= 165 * 0.1
  print(paste0(ct, ':', dim(data[, keep])[2]))
  data = data[, keep]
  data = log2(data + 1)
}

## --------------------- Process clinical metadata table ---------------------
cli = read.csv('/home/lifesci/gchuang/SY/T1D/eQTL/cli/all_cli_v3.csv')
rownames(cli) <- cli$RNA_id
cli = cli[!is.na(cli$geno_map_id),]
bb = data.frame('batch' = unique(cli$batch), 'batch2' = seq(1, length(unique(cli$batch))))
cli = left_join(cli, bb, by = "batch")
rownames(cli) <- cli$RNA_id
cli$sex = as.numeric(cli$sex)
cli$age = as.numeric(cli$age)
cli$batch2 = as.numeric(cli$batch2)

## --------------------- Load cleaned gene location reference ---------------------
ref_clean = read.csv('/home/lifesci/gchuang/SY/T1D/eQTL/gene_ref/gene_ref_clean.csv', row.names = 1)

## --------------------- Load genotype matrix ---------------------
geno = data.frame(fread('./geno/snp_qc_all.csv', header = TRUE))
rownames(geno) = geno$V1
geno = geno[, 2:ncol(geno)]
colnames(geno) = sub(".", "", colnames(geno))

for (ct in c('Pan','CD8T', 'NK', 'CD4T', 'CD14Mono', 'B', 'CD16Mono', 'pDC', 'cDC','T_Prolif', 'Plasmablast')){

  ## --------------------- Gene expression preprocessing ---------------------
  ## Remove genes with zero expression in >90% of individuals, then apply log2 transform
  data = data.frame(fread(paste0('./input/sc-scale/01raw/', ct, '.csv'), header = FALSE))
  var_names = data.frame(t(data[1,]))$X[2:length(data.frame(t(data[1,]))$X)]
  data = data.frame(fread(paste0('./input/sc-scale/01raw/', ct, '.csv')))
  rownames(data) <- data$V1
  data = data[, 2:ncol(data)]
  colnames(data) = var_names

  keep = colSums(data > 0) >= nrow(data) * 0.1
  print(paste0(ct, ':', dim(data[, keep])[2]))
  data = data[, keep]
  data = log2(data + 1)

  ## --------------------- Cell-type-specific gene location info ---------------------
  gene_loc = ref_clean[colnames(data),]

  ## --------------------- Confirm and keep sample IDs for each group ---------------------
  # IDs for T1D and HC
  T1D = as.character(rownames(cli[cli$Disease == 'T1D',]))
  HC  = as.character(rownames(cli[cli$Disease == 'HC',]))

  # Overlap with expression matrix
  T1D_re = intersect(T1D, rownames(data))
  HC_re  = intersect(HC, rownames(data))

  ## --------------------- Subset genotype matrix by group ---------------------
  geno_T1D = geno[, T1D_re]
  geno_HC  = geno[, HC_re]

  ## --------------------- Subset gene expression by group ---------------------
  data_T1D = data.frame(t(data[T1D_re, ]))
  data_HC  = data.frame(t(data[HC_re, ]))
  colnames(data_T1D) <- T1D_re
  colnames(data_HC)  <- HC_re

  ## --------------------- PCA per group ---------------------
  pca_T1D <- data.frame(t(prcomp(t(data_T1D), scale. = FALSE)$x))
  pca_HC  <- data.frame(t(prcomp(t(data_HC),  scale. = FALSE)$x))
  colnames(pca_T1D) <- T1D_re
  colnames(pca_HC)  <- HC_re

  ## --------------------- Covariates per group ---------------------
  cli_T1D = cli[T1D_re, c('sex','age','batch2','Duration')]
  cli_HC  = cli[HC_re,  c('sex','age','batch2')]
  cli_T1D = data.frame(t(cli_T1D))
  cli_HC  = data.frame(t(cli_HC))
  colnames(cli_T1D) = T1D_re
  colnames(cli_HC)  = HC_re

  ## --------------------- Save files per group ---------------------
  ## Check whether output folders exist
  folder_path <- paste0("./input/sc-scale/02QC_log/T1D/", ct)
  if (!file.exists(folder_path)) {
    dir.create(folder_path)
    print(paste("Folder", folder_path, "created successfully."))
  } else {
    print(paste("Folder", folder_path, "already exists."))
  }

  folder_path <- paste0("./input/sc-scale/02QC_log/HC/", ct)
  if (!file.exists(folder_path)) {
    dir.create(folder_path)
    print(paste("Folder", folder_path, "created successfully."))
  } else {
    print(paste("Folder", folder_path, "already exists."))
  }

  ## Save T1D-related outputs
  folder_path <- paste0("./input/sc-scale/02QC_log/T1D/", ct)
  write.table(data_T1D, paste0(folder_path, '/GE.txt'), quote = FALSE)
  write.table(geno_T1D, paste0(folder_path, '/SNP.txt'), quote = FALSE)
  write.table(snploc,   paste0(folder_path, '/snpsloc.txt'), quote = FALSE)
  write.table(cli_T1D,  paste0(folder_path, '/Covariates_cli.txt'), quote = FALSE)
  write.table(pca_T1D,  paste0(folder_path, '/Covariates_pca.txt'), quote = FALSE)
  write.table(gene_loc, paste0(folder_path, '/geneloc.txt'), quote = FALSE)

  ## Save HC-related outputs
  folder_path <- paste0("./input/sc-scale/02QC_log/HC/", ct)
  write.table(data_HC,  paste0(folder_path, '/GE.txt'), quote = FALSE)
  write.table(geno_HC,  paste0(folder_path, '/SNP.txt'), quote = FALSE)
  write.table(snploc,   paste0(folder_path, '/snpsloc.txt'), quote = FALSE)
  write.table(cli_HC,   paste0(folder_path, '/Covariates_cli.txt'), quote = FALSE)
  write.table(pca_HC,   paste0(folder_path, '/Covariates_pca.txt'), quote = FALSE)
  write.table(gene_loc, paste0(folder_path, '/geneloc.txt'), quote = FALSE)
}