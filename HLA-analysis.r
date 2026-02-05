library(dplyr)
library(stringr)
library(ggplot2)
library(vegan)     
library(ggpubr)

clean_hla <- function(x) {
  x <- as.character(x)
  x <- ifelse(substr(x, 1, 1) == "0", substring(x, 2), x)

  parts <- strsplit(x, ":", fixed = TRUE)[[1]]
  if (length(parts) > 2) {
    paste0(parts[1], ":", parts[2])
  } else {
    x
  }
}

# Vectorized wrapper
clean_hla_vec <- function(x) {
  vapply(x, clean_hla, FUN.VALUE = character(1))
}

# -----------------------------
# Helper: compute Shannon diversity for one locus (using allele1 + allele2)
# Returns Shannon index over allele counts.
# -----------------------------
compute_shannon_for_locus <- function(df, locus) {
  a1 <- paste0(locus, ".Allele1")
  a2 <- paste0(locus, ".Allele2")

  # Extract and clean alleles for allele1
  test1 <- df[, c(a1, "disease_status", "sex")]
  colnames(test1) <- c("HLA", "disease_status", "sex")
  test1$HLA <- clean_hla_vec(test1$HLA)

  # Extract and clean alleles for allele2
  test2 <- df[, c(a2, "disease_status", "sex")]
  colnames(test2) <- c("HLA", "disease_status", "sex")
  test2$HLA <- clean_hla_vec(test2$HLA)

  # Dummy encode alleles (one-hot per allele) for each allele column
  dummy1 <- as.data.frame(model.matrix(~0 + HLA, data = test1))
  dummy2 <- as.data.frame(model.matrix(~0 + HLA, data = test2))

  # Align columns and sum counts across allele1+allele2
  common_cols <- intersect(colnames(dummy1), colnames(dummy2))
  inte <- dummy1[, common_cols, drop = FALSE] + dummy2[, common_cols, drop = FALSE]

  only1 <- setdiff(colnames(dummy1), common_cols)
  only2 <- setdiff(colnames(dummy2), common_cols)

  # Combine all allele counts into one matrix
  if (length(only1) > 0 && length(only2) > 0) {
    mat <- cbind(inte, dummy1[, only1, drop = FALSE], dummy2[, only2, drop = FALSE])
  } else if (length(only1) > 0 && length(only2) == 0) {
    mat <- cbind(inte, dummy1[, only1, drop = FALSE])
  } else if (length(only1) == 0 && length(only2) > 0) {
    mat <- cbind(inte, dummy2[, only2, drop = FALSE])
  } else {
    mat <- inte
  }

  # Allele counts across all individuals
  allele_counts <- colSums(mat)

  # Shannon diversity index
  shannon_index <- vegan::diversity(allele_counts, index = "shannon")

  data.frame(HLA = locus, shannon_index = shannon_index, stringsAsFactors = FALSE)
}

# -----------------------------
# Compute Shannon indices for a cohort (T1D or HC)
# -----------------------------
compute_shannon_df <- function(file) {
  dat <- read.csv(file, row.names = 1, stringsAsFactors = FALSE)

  loci <- c("A","B","C","DRB1","DQA1","DQB1","DPA1","DPB1")
  out <- do.call(rbind, lapply(loci, function(l) compute_shannon_for_locus(dat, l)))
  out
}

# -----------------------------
# Run for T1D and HC
# -----------------------------
df_T1D <- compute_shannon_df("./T1D-HLA-allele-with-cli.csv")
df_T1D$disease <- "T1D"

df_HC <- compute_shannon_df("./HC-HLA-allele-with-cli.csv")
df_HC$disease <- "HC"

ALL_shannon <- rbind(df_T1D, df_HC)
ALL_shannon$disease <- factor(ALL_shannon$disease, levels = c("HC", "T1D"))

# -----------------------------
# Plot: grouped boxplot + dotplot + paired test
# Note: color_m must be defined elsewhere (a named vector for HLA fill colors).
# -----------------------------
p3 <- ggplot(ALL_shannon, aes(x = disease, y = shannon_index)) +
  geom_boxplot() +
  theme_bw() +
  geom_dotplot(
    aes(fill = HLA),
    trim = FALSE,
    binaxis = "y",
    stackdir = "centerwhole",
    dotsize = 1.2,
    position = position_dodge(0.2)
  ) +
  scale_fill_manual(values = color_m) +
  labs(
    title = "",
    x = "Number of digits of HLA alleles",
    y = "Shannon's Diversity Index"
  ) +
  geom_line(aes(group = HLA), color = "gray", linewidth = 0.5) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.4),
    axis.text.y = element_text(color = "black", size = 10),
    axis.text.x = element_text(margin = margin(t = 2), color = "black", size = 10),
    legend.position = "none",
    panel.spacing = unit(0, "lines"),
    panel.grid = element_blank()
  ) +
  ggpubr::stat_compare_means(
    method = "t.test",
    paired = TRUE,
    comparisons = list(c("T1D", "HC"))
  ) +
  ylim(0.5, 3.6)

p3
                               
# Logistic regression
# Keep alleles with MAF > 0.01 (implemented as count >= 1% of samples later)
merge_df2 = rbind(T1D,HC)                              
HLA_res_all <- data.frame()

for (i in c('A','B','C','DRB1','DQA1','DQB1','DPA1','DPB1')) {

  i1 <- paste0(i, '.Allele1')
  i2 <- paste0(i, '.Allele2')

  # -----------------------------
  # Allele1: clean allele strings and one-hot encode
  # -----------------------------
  test <- merge_df2[c(i1, 'disease_status', 'sex', 'age', 'province2')]
  colnames(test) <- c('HLA', 'disease_status', 'sex', 'age', 'province2')

  vec <- test$HLA
  vec <- ifelse(substr(vec, 1, 1) == "0", substring(vec, 2), vec)
  test$HLA <- vec

  # Reduce allele resolution to 2-field
  re <- c()
  for (x in seq_along(test$HLA)) {
    parts <- strsplit(as.character(test$HLA[x]), ":", fixed = TRUE)[[1]]
    re_x <- ifelse(length(parts) > 2, paste0(parts[1], ":", parts[2]), test$HLA[x])
    re <- append(re, re_x)
  }
  test$HLA <- re

  dummy_data <- model.matrix(~ 0 + HLA, data = test)
  dummy_data <- data.frame(dummy_data)

  # -----------------------------
  # Allele2: clean allele strings and one-hot encode
  # -----------------------------
  test <- merge_df2[c(i2, 'disease_status', 'sex', 'age', 'province2')]
  colnames(test) <- c('HLA', 'disease_status', 'sex', 'age', 'province2')

  vec <- test$HLA
  vec <- ifelse(substr(vec, 1, 1) == "0", substring(vec, 2), vec)
  test$HLA <- vec

  re <- c()
  for (x in seq_along(test$HLA)) {
    parts <- strsplit(as.character(test$HLA[x]), ":", fixed = TRUE)[[1]]
    re_x <- ifelse(length(parts) > 2, paste0(parts[1], ":", parts[2]), test$HLA[x])
    re <- append(re, re_x)
  }
  test$HLA <- re

  dummy_data2 <- model.matrix(~ 0 + HLA, data = test)
  dummy_data2 <- data.frame(dummy_data2)

  # -----------------------------
  # Combine allele1 + allele2 encodings
  # -----------------------------
  v1 <- intersect(colnames(dummy_data), colnames(dummy_data2))
  v2 <- colnames(dummy_data[, colSums(dummy_data) >= 3, drop = FALSE])
  v3 <- colnames(dummy_data2[, colSums(dummy_data2) >= 3, drop = FALSE])

  inte <- dummy_data[, v1, drop = FALSE] + dummy_data2[, v1, drop = FALSE]
  a1 <- dummy_data[, setdiff(v2, v1), drop = FALSE]
  a2 <- dummy_data2[, setdiff(v3, v1), drop = FALSE]

  if (ncol(a1) > 0 && ncol(a2) > 0) {
    test <- cbind(inte, a1, a2)
  } else if (ncol(a1) > 0 && ncol(a2) == 0) {
    test <- cbind(inte, a1)
  } else if (ncol(a1) == 0 && ncol(a2) > 0) {
    test <- cbind(inte, a2)
  } else {
    test <- inte
  }

  # Filter rare alleles: keep columns with count >= 1% of total samples
  test <- test[, colSums(test) >= nrow(merge_df2) * 0.01, drop = FALSE]

  test$sex <- as.numeric(merge_df2$sex)
  test$disease_status <- as.numeric(merge_df2$disease_status)
  test$age <- as.numeric(merge_df2$age)
  test$province2 <- as.numeric(merge_df2$province2)

  # -----------------------------
  # Fit logistic regression per allele
  # -----------------------------
  name <- c()
  p <- c()
  coe <- c()
  OR <- c()

  # Iterate over allele columns (excluding covariates + PCs)
  n_cov <- length(c(colnames(pc), 'sex', 'disease_status', 'age', 'province2'))
  for (j in seq_len(ncol(test) - n_cov)) {

    temp <- test[, c(colnames(test)[j], 'sex', 'disease_status', 'age', 'province2', colnames(pc)), drop = FALSE]
    colnames(temp) <- c('loci', 'sex', 'disease_status', 'age', 'province2', colnames(pc))

    model <- glm(disease_status ~ loci + sex + province2, family = binomial(link = "logit"), data = temp)
    null_model <- glm(disease_status ~ sex + province2, family = binomial(link = "logit"), data = temp)

    likelihood_ratio_test <- anova(null_model, model, test = "Chisq")
    p_value <- likelihood_ratio_test$Pr[2]

    or <- exp(coef(model))['loci']

    name <- append(name, colnames(test)[j])
    p <- append(p, p_value)
    OR <- append(OR, or)
    coe <- append(coe, coef(model)["loci"])
  }


  HLA_res <- data.frame(name, p, coe, OR)
  HLA_res$loci <- i
  HLA_res_all <- rbind(HLA_res_all, HLA_res)
}

                               
## Linear regression with age at diagnosis (onset_age)
cli <- read.csv('T1D_cli.csv')
merge_df3 <- merge(T1D, cli[c('sample_id','onset_age')], by = 'sample_id')

HLA_res_all <- data.frame()

for (i in c('A','B','C','DRB1','DQA1','DQB1','DPA1','DPB1')) {

  i1 <- paste0(i, '.Allele1')
  i2 <- paste0(i, '.Allele2')

  # -----------------------------
  # Allele1: clean allele strings and one-hot encode
  # -----------------------------
  test <- merge_df3[c(i1, 'disease_status', 'sex', 'age', 'province2')]
  colnames(test) <- c('HLA', 'disease_status', 'sex', 'age', 'province2')

  # Remove leading '0' in allele codes
  vec <- test$HLA
  vec <- ifelse(substr(vec, 1, 1) == "0", substring(vec, 2), vec)
  test$HLA <- vec

  # Reduce allele resolution to 2-field (e.g., 01:02:03 -> 01:02)
  re <- c()
  for (x in seq_along(test$HLA)) {
    parts <- strsplit(as.character(test$HLA[x]), ":", fixed = TRUE)[[1]]
    re_x <- ifelse(length(parts) > 2, paste0(parts[1], ":", parts[2]), test$HLA[x])
    re <- append(re, re_x)
  }
  test$HLA <- re

  dummy_data <- model.matrix(~ 0 + HLA, data = test)
  dummy_data <- data.frame(dummy_data)

  # -----------------------------
  # Allele2: clean allele strings and one-hot encode
  # -----------------------------
  test <- merge_df3[c(i2, 'disease_status', 'sex', 'province2')]
  colnames(test) <- c('HLA', 'disease_status', 'sex', 'province2')

  vec <- test$HLA
  vec <- ifelse(substr(vec, 1, 1) == "0", substring(vec, 2), vec)
  test$HLA <- vec

  re <- c()
  for (x in seq_along(test$HLA)) {
    parts <- strsplit(as.character(test$HLA[x]), ":", fixed = TRUE)[[1]]
    re_x <- ifelse(length(parts) > 2, paste0(parts[1], ":", parts[2]), test$HLA[x])
    re <- append(re, re_x)
  }
  test$HLA <- re

  dummy_data2 <- model.matrix(~ 0 + HLA, data = test)
  dummy_data2 <- data.frame(dummy_data2)

  # -----------------------------
  # Combine allele1 + allele2 encodings
  # -----------------------------
  v1 <- intersect(colnames(dummy_data), colnames(dummy_data2))
  v2 <- colnames(dummy_data[, colSums(dummy_data) >= 3, drop = FALSE])
  v3 <- colnames(dummy_data2[, colSums(dummy_data2) >= 3, drop = FALSE])

  inte <- dummy_data[, v1, drop = FALSE] + dummy_data2[, v1, drop = FALSE]
  a1 <- dummy_data[, setdiff(v2, v1), drop = FALSE]
  a2 <- dummy_data2[, setdiff(v3, v1), drop = FALSE]

  if (ncol(a1) > 0 && ncol(a2) > 0) {
    test <- cbind(inte, a1, a2)
  } else if (ncol(a1) > 0 && ncol(a2) == 0) {
    test <- cbind(inte, a1)
  } else if (ncol(a1) == 0 && ncol(a2) > 0) {
    test <- cbind(inte, a2)
  } else {
    test <- inte
  }

  # Filter rare alleles: keep columns with count >= 1% of total samples
  test <- test[, colSums(test) >= nrow(merge_df3) * 0.01, drop = FALSE]

  # -----------------------------
  # Add covariates for regression
  # -----------------------------
  test$onset_age <- merge_df3$onset_age
  test$sex <- as.numeric(merge_df3$sex)
  test$age <- as.numeric(merge_df3$age)
  test$province2 <- as.numeric(merge_df3$province2)

  # -----------------------------
  # Fit linear regression per allele
  # -----------------------------
  name <- c()
  p <- c()
  coe <- c()

  # Loop over allele columns (exclude onset_age/sex/age/province2 => 4 columns)
  for (j in seq_len(ncol(test) - 4)) {

    temp <- test[, c(colnames(test)[j], 'sex', 'onset_age', 'province2'), drop = FALSE]
    colnames(temp) <- c('loci', 'sex', 'onset_age', 'province2')

    model <- lm(onset_age ~ loci + sex + province2, data = temp)
    null_model <- lm(onset_age ~ sex + province2, data = temp)

    F <- anova(null_model, model)
    p_value <- F$`Pr(>F)`[2]

    name <- append(name, colnames(test)[j])
    p <- append(p, p_value)
    coe <- append(coe, coef(model)["loci"])
  }

  HLA_res <- data.frame(name = name, p = p, coe = coe)
  HLA_res$loci <- i
  HLA_res_all <- rbind(HLA_res_all, HLA_res)
}
