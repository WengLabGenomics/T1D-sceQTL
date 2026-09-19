library(data.table)
library(dplyr)
library(coloc)

cts = c('CD14Mono', 'CD16Mono', 'cDC', 'pDC','CD4T','CD8T', 'NK','T_Prolif', 'B', 'Plasma')
for (ct in cts) {
geno = data.frame(fread(paste0('./input/disease_all_geno/',ct,'_GE.txt')))
rownames(geno) =geno$V1
geno = geno[,2:ncol(geno)]
colnames(geno) = sub(".", "", colnames(geno))
    sdY_all <- apply(geno , 1, sd, na.rm = TRUE)
    sdY_df <- data.frame(
      gene = names(sdY_all),
      sdY = as.numeric(sdY_all)
    )
    write.csv(sdY_df,paste0('./input/disease_all_geno/',ct,'_sdY.csv'))
    }

risk <- read.csv("./res/T1D_risk_eQTL_all.csv")

for (ct in cts) {
    x = readRDS(paste0('./res/combine-All/',ct,'_eQTL.rds'))
    x2 = me$cis$eqtls
    x3 = x2[x2$gene %in% risk$gene,]
    write.csv(x3, paste0("./res/risk_eGene/", ct, "_cis_all.csv"))
}

gwas <- data.frame(fread('./GWAS_sum/2021-NG-hg38.ma'))
n_case = 16159
n_control = 25386
n_eqtl = 163

all_result <- list()
for (ct in cts) {

    message("Processing: ", ct)

    re_g <- unique(risk$gene[risk$ct == ct])

    eqtl_0 <- fread(
        paste0("./res/risk_eGene/", ct, "_cis_all.csv")
    )

    eqtl_0$se <- abs(eqtl_0$beta / eqtl_0$statistic)

    sdY_df <- read.csv(
        paste0("./input/disease_all_geno/", ct, "_sdY.csv")
    )
    rownames(sdY_df) <- sdY_df$gene

    result_list <- list()
    snp_result_list <- list()

    for (gene_i in re_g) {

        eqtl <- eqtl_0[gene == gene_i]

        if (nrow(eqtl) == 0) next


        res <- tryCatch(

            coloc.abf(

                dataset1 = list(
                    beta = gwas$b,
                    varbeta = gwas$se^2,
                    snp = gwas$SNP,
                    type = "cc",
                    s = n_case / (n_case + n_control),
                    N = n_case + n_control
                ),

                dataset2 = list(
                    beta = eqtl$beta,
                    varbeta = eqtl$se^2,
                    snp = eqtl$snps,
                    sdY = sdY_df[gene_i, "sdY"],
                    type = "quant",
                    N = n_eqtl
                )

            ),

            error = function(e) {

                message(
                    "Skip ", ct, " | ", gene_i,
                    " | ", e$message
                )

                return(NULL)
            }
        )


        if (is.null(res)) next


        snp_res <- as.data.frame(res$results)

        top_idx <- which.max(snp_res$SNP.PP.H4)

        top_shared_SNP <- snp_res$snp[top_idx]
        top_SNP_PP_H4 <- snp_res$SNP.PP.H4[top_idx]


        PP3 <- as.numeric(res$summary["PP.H3.abf"])
        PP4 <- as.numeric(res$summary["PP.H4.abf"])

        result_list[[gene_i]] <- data.frame(

            ct = ct,
            gene = gene_i,

            nsnp = as.numeric(
                res$summary["nsnps"]
            ),

            PP0 = as.numeric(
                res$summary["PP.H0.abf"]
            ),

            PP1 = as.numeric(
                res$summary["PP.H1.abf"]
            ),

            PP2 = as.numeric(
                res$summary["PP.H2.abf"]
            ),

            PP3 = PP3,
            PP4 = PP4,

            PP4_conditional = PP4 / (PP3 + PP4),

            top_shared_SNP = top_shared_SNP,

            top_SNP_PP_H4 = top_SNP_PP_H4,

            status = "success"
        )

        snp_res$ct <- ct
        snp_res$gene <- gene_i

        snp_result_list[[gene_i]] <- snp_res

    }

    coloc_result <- bind_rows(result_list)

    coloc_result <- coloc_result %>%
        arrange(desc(PP4))


    fwrite(
        coloc_result,
        paste0(
            "./coloc/",
            ct,
            "_coloc_summary.csv"
        )
    )

    coloc_snp_result <- bind_rows(
        snp_result_list
    )

    fwrite(
        coloc_snp_result,
        paste0(
            "./coloc/2021-NG/",
            ct,
            "_coloc_SNP.csv.gz"
        )
    )

    all_result[[ct]] <- coloc_result

}

all_coloc_result <- bind_rows(
    all_result
) %>%
    arrange(desc(PP4))


fwrite(
    all_coloc_result,
    "./coloc/2021-NG/all_celltypes_coloc_summary.csv"
)