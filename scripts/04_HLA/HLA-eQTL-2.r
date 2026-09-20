library(ggplot2)
library(dplyr)
library(RColorBrewer)
library("data.table")
library(stringr)

hr = read.csv('./HLA_res/disease_sig_hla_loci.csv')
hr = hr[hr$p < 5e-8,]
df_all = data.frame()
for (ct in c('CD14Mono', 'CD16Mono', 'cDC', 'pDC','CD4T','CD8T', 'NK','T_Prolif', 'B', 'Plasma')) {
    df =  data.frame(fread(paste0("./res/HLA_classic_All/",ct,"_cis_fdr005.csv")))
    df$ct = ct
    df_all = rbind(df_all,df)
    }

df_all = df_all[df_all$snps %in% hr$name,]

library(stringr)
x1 = str_split_fixed(df_all$snps,'_',3)
x2 = str_split_fixed(x1[,3],'\\.',2)
re = sapply(x2[,1], function(x) {
 if (grepl("^\\d{2,}$", x)) {
    return(x)  # 如果已经是两个数字，保持不变
  } else {
    return(paste0("0", x))  # 如果不是两个数字，补充一个0
  }
})
names(re) <- NULL
df_all$snps = paste0(x1[,2],'*',re,':',x2[,2])

df_all = df_all[c('ct','gene','snps','beta','statistic','pvalue','FDR')]

df_all$eqtl = paste0(df_all$snps,'_',df_all$gene)

simes_test <- function(pvals){
  
  # remove NA
  pvals <- pvals[!is.na(pvals)]
  
  if(length(pvals) == 0){
    return(NA)
  }
  
  # sort p-values
  pvals <- sort(pvals)
  
  m <- length(pvals)
  
  # Simes
  simes_p <- min(pvals * m / seq_len(m))
  
  return(min(simes_p, 1))
}

library(dplyr)
simes_result <- df_all %>%
  filter(!is.na(snps),
         !is.na(ct),
         !is.na(pvalue)) %>%
  group_by(snps, ct) %>%
  summarise(
    n_gene = n(),
    simes_p = simes_test(pvalue),
    .groups = "drop"
  )

simes_result <- simes_result %>%
  mutate(
    simes_FDR = p.adjust(
      simes_p,
      method = "BH"
    )
  )

library(dplyr)
library(tidyr)
library(pheatmap)
library(tibble)

heatmap_df <- simes_result %>%
  filter(!is.na(snps),
         !is.na(ct),
         !is.na(simes_FDR)) %>%
  mutate(
    score = -log10(simes_FDR)
  ) %>%
  select(snps, ct, score)


# 转换为matrix
heatmap_mat <- heatmap_df %>%
  pivot_wider(
    names_from = ct,
    values_from = score
  ) %>%
  column_to_rownames("snps") %>%
  as.matrix()
heatmap_mat[is.na(heatmap_mat)] =0
heatmap_mat[is.infinite(heatmap_mat)] <- max(
  heatmap_mat[is.finite(heatmap_mat)]
)
heatmap_mat[heatmap_mat > 15] <- 15

ct_order  = c('CD14Mono', 'CD16Mono', 'cDC', 'pDC','CD4T','CD8T', 'NK','T_Prolif', 'B', 'Plasma')
allele_order <- c(
  grep("^A\\*", rownames(heatmap_mat), value=TRUE),
  grep("^B\\*", rownames(heatmap_mat), value=TRUE),
  grep("^C\\*", rownames(heatmap_mat), value=TRUE),
  
  grep("^DR", rownames(heatmap_mat), value=TRUE),
  grep("^DQ", rownames(heatmap_mat), value=TRUE),
  grep("^DP", rownames(heatmap_mat), value=TRUE)
)

heatmap_mat <- heatmap_mat[allele_order,]

heatmap_mat <- heatmap_mat[,ct_order]

row_annotation <- data.frame(
  Class = ifelse(
    grepl("^(A|B|C)", rownames(heatmap_mat)),
    "Class I",
    "Class II"
  )
)

rownames(row_annotation) <- rownames(heatmap_mat)

annotation_colors <- list(
  Class = c(
    "Class I" = "#9ECAE1",
    "Class II" = "#BC80BD"
  )
)

pheatmap(
  t(heatmap_mat),
  annotation_col = row_annotation,
      annotation_colors = annotation_colors,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("white", "#FEE8C8", "#FDBB84", "#E34A33", "#B30000"))(100),
  border_color = 'white',
    cellwidth = 15,
  cellheight = 15)

df_all$class = 'Class-II'
df_all[grepl("^A\\*",df_all$snps) | grepl("^B\\*",df_all$snps)|  grepl("^C\\*",df_all$snps),'class'] = 'Class-I'

simes_result <- df_all %>%
  filter(!is.na(class),
         !is.na(ct),
         !is.na(pvalue)) %>%
  group_by(class, ct) %>%
  summarise(
    n_gene = n(),
    simes_p = simes_test(pvalue),
    .groups = "drop"
  )

simes_result <- simes_result %>%
  mutate(
    simes_FDR = p.adjust(
      simes_p,
      method = "BH"
    )
  )

# 计算heatmap score
heatmap_df <- simes_result %>%
  filter(!is.na(class),
         !is.na(ct),
         !is.na(simes_FDR)) %>%
  mutate(
    score = -log10(simes_FDR)
  ) %>%
  select(class, ct, score)


# 转换为matrix
heatmap_mat <- heatmap_df %>%
  pivot_wider(
    names_from = ct,
    values_from = score
  ) %>%
  column_to_rownames("class") %>%
  as.matrix()
heatmap_mat[is.na(heatmap_mat)] =0
heatmap_mat[is.infinite(heatmap_mat)] <- max(
  heatmap_mat[is.finite(heatmap_mat)]
)
heatmap_mat[heatmap_mat > 15] <- 15
heatmap_mat <- heatmap_mat[,ct_order]

pheatmap(
  t(heatmap_mat),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("white", "#FEE8C8", "#FDBB84", "#E34A33", "#B30000"))(100),
  border_color = 'white',
    cellwidth = 15,
  cellheight = 15)