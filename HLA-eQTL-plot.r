library(data.table)
library(stringr)
library(dplyr)

library(ggplot2)
library(RColorBrewer)
library(ggrastr) 
library(ggrepel)
library(patchwork) 

# -----------------------------
# Inputs
# -----------------------------
loci <- c('A','B','C','DRB1','DQA1','DQB1','DPB1','DPA1')
pos  <- c(29944214,31355516,31270415,32584287,32640529,32663284,33081591,33069309)
df_loci <- data.frame(loci = loci, pos = pos, stringsAsFactors = FALSE)

bim <- read.table('./chr6_imputed/chr6_imputed.bim', stringsAsFactors = FALSE)
bim$snps <- gsub('-', '.', paste0(bim$V2, '_', bim$V5))
bim$pos  <- bim$V4

gtf_hla <- read.table('./gene_ref/hla.txt', header = TRUE, stringsAsFactors = FALSE)

gene_color <- brewer.pal(8, 'Set1')

dis <- 'T1D'
celltypes <- c('Pan','CD14Mono','CD16Mono','cDC','pDC','CD4T','CD8T','NK','T_Prolif','B','Plasmablast')

keep_genes <- c('HLA-A','HLA-B','HLA-C','HLA-DQA1','HLA-DQB1','HLA-DRB1','HLA-DPA1','HLA-DPB1')

xlim_min <- 28.002253
xlim_max <- 33.985869

for (ct in celltypes) {

  # -----------------------------
  # Read results
  # -----------------------------
  results1 <- read.csv(paste0('./res/HLA-imputed/', dis, '/', ct, '_fdr005.csv'), row.names = 1, stringsAsFactors = FALSE)
  results2 <- read.csv(paste0('./res/HLA-classical/', dis, '/', ct, '_cis_fdr005.csv'), row.names = 1, stringsAsFactors = FALSE)

  # Add position for imputed SNPs
  results1 <- left_join(results1, bim[, c('snps','pos')], by = 'snps')

  # Add position for classical alleles (from locus mapping)
  # assumes results2$snps looks like "..._<LOCUS>_..." and locus is in the 2nd field
  locus2 <- str_split_fixed(results2$snps, '_', 3)[, 2]
  results2$pos <- df_loci$pos[match(locus2, df_loci$loci)]

  # Merge
  results <- bind_rows(results1, results2) %>%
    mutate(
      pos = as.numeric(pos),
      cell_type = ct
    ) %>%
    arrange(pos) %>%
    filter(gene %in% keep_genes)

  # Lead variant per (cell_type, gene)
  lead_variants <- results %>%
    group_by(cell_type, gene) %>%
    slice_min(order_by = pvalue, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(variant = snps)

  # -----------------------------
  # Plot
  # -----------------------------
  p1 <- results %>%
    sample_frac(1) %>%  # shuffle to reduce overplotting bias
    ggplot() +
    geom_point_rast(aes(x = pos / 1e6, y = -log10(pvalue), col = gene), size = 0.3) +
    ylab('-log10(P-value)') +
    xlab('Chr 6 position (MB)') +
    theme_classic() +
    scale_color_manual(values = gene_color) +
    theme(legend.position = 'none') +
    guides(colour = guide_legend(override.aes = list(size = 4))) +
    ggrepel::geom_text_repel(
      data = lead_variants,
      aes(x = pos / 1e6, y = -log10(pvalue), label = variant),
      size = 3,
      segment.size = 0.3,
      min.segment.length = 0,
      nudge_y = 10,
      nudge_x = 0.1,
      max.overlaps = Inf
    ) +
    xlim(xlim_min, xlim_max)

  p2 <- ggplot(gtf_hla) +
    geom_point(aes(x = TSS / 1e6, y = plot_pos, col = hla_gene), size = 3, shape = 18) +
    scale_color_manual(values = gene_color) +
    theme_void() +
    theme(legend.position = 'none') +
    xlim(xlim_min, xlim_max) +
    ylim(-1, 1)

  print(ct)
  print(p1 / p2 + plot_layout(heights = c(5, 1)))
}