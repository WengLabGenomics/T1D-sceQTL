library(MatrixEQTL)

library(MatrixEQTL)

# -----------------------------
# Script arguments:
#   1) Disease group: T1D / HC
#   2) Cell type
#   3) Number of genotype PCs to include
# -----------------------------
arg <- commandArgs(TRUE)
dis <- arg[1]
ct  <- arg[2]
n   <- as.numeric(arg[3])

# -----------------------------
# Paths
# -----------------------------
base.dir <- './demo/'
save.dir <- './res/'

# Covariates to include (clinical)
re_cov_T1D <- c('sex','age','batch2','Duration')
re_cov_HC  <- c('sex','age','batch2')

c1 = read.table(paste0(base.dir,dis,'/',ct, "/Covariates_cli.txt"))
# Genotype PCs (rows = PCs, cols = samples)
c3 <- read.table(paste0(base.dir, dis, "/", dis, "_geno_pca.txt"),
                 header = TRUE, check.names = FALSE)

# Select genotype PCs to include
cov_pc <- rownames(c3)[1:n]

# Build covariate row list
if (dis == 'T1D') {
  re_cov <- append(re_cov_T1D, cov_pc)
} else {
  re_cov <- append(re_cov_HC, cov_pc)
}

cov = rbind(c1,c3[colnames(c1)])
cov = cov[re_cov,]
write.table(cov,paste0(base.dir,dis,'/',ct, "/Covariates.txt"),quote = FALSE) q

# -----------------------------
# MatrixEQTL settings
# -----------------------------
useModel <- modelLINEAR

SNP_file_name             <- paste0(base.dir, dis, '/', ct, "/SNP.txt")
snps_location_file_name   <- paste0(base.dir, dis, '/', ct, "/snpsloc.txt")

expression_file_name      <- paste0(base.dir, dis, '/', ct, "/GE.txt")
gene_location_file_name   <- paste0(base.dir, dis, '/', ct, "/geneloc.txt")

covariates_file_name      <- paste0(base.dir, dis, '/', ct, "/Covariates.txt")

output_file_name_cis <- tempfile()
output_file_name_tra <- tempfile()

pvOutputThreshold_cis <- 2e-2
pvOutputThreshold_tra <- 1e-2

errorCovariance <- numeric()
cisDist <- 1e6

# -----------------------------
# Load genotype data
# -----------------------------
snps <- SlicedData$new()
snps$fileDelimiter       <- " "
snps$fileOmitCharacters   <- "NA"
snps$fileSkipRows         <- 1
snps$fileSkipColumns      <- 1
snps$fileSliceSize        <- 2000
snps$LoadFile(SNP_file_name)

# -----------------------------
# Load gene expression data
# -----------------------------
gene <- SlicedData$new()
gene$fileDelimiter       <- " "
gene$fileOmitCharacters   <- "NA"
gene$fileSkipRows         <- 1
gene$fileSkipColumns      <- 1
gene$fileSliceSize        <- 2000
gene$LoadFile(expression_file_name)

# -----------------------------
# Load covariates
# -----------------------------
cvrt <- SlicedData$new()
cvrt$fileDelimiter       <- " "
cvrt$fileOmitCharacters   <- "NA"
cvrt$fileSkipRows         <- 1
cvrt$fileSkipColumns      <- 1
if (length(covariates_file_name) > 0) {
  cvrt$LoadFile(covariates_file_name)
}

# -----------------------------
# SNP / gene positions
# -----------------------------
snpspos <- read.table(snps_location_file_name,
                      header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
snpspos$snpid <- rownames(snpspos)
snpspos <- snpspos[, c('snpid','chr','pos')]

genepos <- read.table(gene_location_file_name,
                      header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
genepos$geneid <- rownames(genepos)
genepos <- genepos[, c('geneid','chr','left','right')]

# -----------------------------
# Run MatrixEQTL
# -----------------------------
me <- Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = output_file_name_tra,
  pvOutputThreshold = pvOutputThreshold_tra,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = TRUE,
  noFDRsaveMemory = FALSE
)

unlink(output_file_name_tra)
unlink(output_file_name_cis)

# -----------------------------
# Save results
# -----------------------------
n_n <- paste0('new-genoPC-', as.character(n))
folder_path <- paste0(save.dir, dis, '/', n_n, '/')
dir.create(folder_path, recursive = TRUE, showWarnings = FALSE)

saveRDS(me, paste0(folder_path, ct, "_eQTL.rds"))

# Keep only cis results with FDR < 0.05
res <- me$cis$eqtls
res <- res[res$FDR < 0.05, ]
write.csv(res, paste0(folder_path, ct, "_cis_eQTL_fdr005.csv"), row.names = FALSE)

# -----------------------------
# Logging
# -----------------------------
print(paste0('Now is: Disease_state-', dis, '; Celltype-', ct))
print(paste0('Analysis done in: ', me$time.in.sec, ' seconds'))
                                                                   
