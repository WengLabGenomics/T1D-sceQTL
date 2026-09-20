library(tidyverse)
library(data.table)
library(lme4)
library(Matrix)
library(edgeR)
library(DESeq2)
library(Hmisc)

cli  = read.csv('./input/all_cli_v3.csv')
rownames(cli) <- cli$RNA_id
for (ct in c('CD8T', 'NK', 'CD4T', 'CD14Mono', 'B', 'CD16Mono', 'pDC', 'cDC','T_Prolif', 'Plasma')){
data = data.frame(fread(paste0('./input/01raw/',ct,'.csv'),header = FALSE))
var_names = data.frame(t(data[1,]))$X[2:length(data.frame(t(data[1,]))$X)]
data = data.frame(fread(paste0('./input/01raw/',ct,'.csv')))
rownames(data) <- data$V1
data = data[,2:ncol(data)]
colnames(data) = var_names
keep = colSums(data > 0 ) >= nrow(data)*0.1
print(paste0(ct,':',dim(data[,keep])[2]))
data = data[,keep]
data  = log2(data+1)
T1D = as.character(rownames(cli[cli$Disease == 'T1D',]))
HC = as.character(rownames(cli[cli$Disease == 'HC',]))
T1D_re = intersect(T1D,rownames(data))
HC_re = intersect(HC,rownames(data))
data_T1D = data.frame(t(data[T1D_re,]))
data_HC = data.frame(t(data[HC_re,]))
colnames(data_T1D) <- T1D_re
colnames(data_HC) <- HC_re
write.csv(data_T1D,paste0('./cli_linear/gene_exp/sc-scale/T1D/',ct,'_exp.csv'),quote = FALSE)
write.csv(data_HC,paste0('./cli_linear/gene_exp/sc-scale/HC/',ct,'_exp.csv'),quote = FALSE)
    }

HC_cli  =  cli[cli$Disease == 'HC',]
T1D_cli  =  cli[cli$Disease == 'T1D',]

write.csv(HC_cli ,'./clinical_data/HC_clinical_data/HC_cli_final.csv',quote = FALSE)
write.csv(T1D_cli ,'./clinical_data/T1D_clinical_data/T1D_cli_final.csv',quote = FALSE)

for (Disease in c('HC','T1D')){
    for (ct in c('CD8T', 'NK', 'CD4T', 'CD14Mono', 'B', 'CD16Mono', 'pDC', 'cDC','T_Prolif', 'Plasma')){
        exp <- data.frame(fread(paste0('./cli_linear/gene_exp/sc-scale/',Disease,'/',clustering,'_exp.csv'),header = TRUE))
        rownames(exp) <- exp$V1
        exp = exp[,2:ncol(exp)]
        colnames(exp) <- substr(colnames(exp),2,nchar(colnames(exp)))

        if (Disease == 'HC') {
            cli <- read.csv('./clinical_data/HC_clinical_data/HC_cli_final.csv',row.names = 1)
            rownames(cli) <- as.character(rownames(cli))
            sample <- intersect(colnames(exp),rownames(cli))
            cli = cli[sample,]
            gene_exp = exp[,sample]
            all_out = data.frame()
            for (gene in rownames(gene_exp)){
                E <- as.numeric(gene_exp[gene,])
                age <- as.numeric(cli$age)
                data <- data.frame(E, age)
                res <- rcorr(as.matrix(data),type = "spearman")
                r = res$r[1,2]
                p = res$P[1,2]
                out = data.frame(gene = gene,r = r,p=p,ct = clustering )
                all_out = rbind(all_out,out)
                }
                write.csv(all_out,paste0('./cli_linear/res/HC_cor/',clustering,'_age_cor.csv'))
        } else if (Disease == 'T1D') {
           cli <- read.csv('./clinical_data/T1D_clinical_data/T1D_cli_final.csv',row.names = 1)
            rownames(cli) <- as.character(rownames(cli))
            sample <- intersect(colnames(exp),rownames(cli))
            cli = cli[sample,]
            gene_exp = exp[,sample]
            all_out = data.frame()
            for (gene in rownames(gene_exp)){
                E <- as.numeric(gene_exp[gene,])
                onset_age <- as.numeric(cli$onset_age)
                data <- data.frame(E,onset_age)
                res <- rcorr(as.matrix(data),type = "spearman")
                r = res$r[1,2]
                p = res$P[1,2]
                out = data.frame(gene = gene,r = r,p=p,ct = clustering )
                all_out = rbind(all_out,out)
                }
            write.csv(all_out,paste0('./cli_linear/res/T1D_cor/',clustering,'_age_cor.csv'))
        } else{
            print("Don not forget the disease status!!")
            }
    }
}

eg_l = c()
all = data.frame()
for (ct in c('CD8T', 'NK', 'CD4T', 'CD14Mono', 'B', 'CD16Mono', 'pDC', 'cDC','T_Prolif', 'Plasma')) {
    dis = 'T1D'
    res = read.csv(paste0('./cli_linear/res/',dis,'_cor/',ct,'_age_cor.csv'))
    res$fdr <- p.adjust(res$p, method = "fdr")
    de <- res[!is.na(res$fdr) & res$fdr < 0.05 & res$r < 0,]$gene
    dl <- res[!is.na(res$fdr) & res$fdr < 0.05 & res$r > 0,]$gene
    dis = 'HC'
    res2 = read.csv(paste0('./cli_linear/res/',dis,'_cor/',ct,'_age_cor.csv'))
    res2 = res2[!is.na(res2$p),]
    res2$fdr <- p.adjust(res2$p, method = "fdr")
    he <- res2[!is.na(res2$fdr) & res2$fdr < 0.05 & res2$r < 0,]$gene
    hl <- res2[!is.na(res2$fdr) & res2$fdr < 0.05 & res2$r > 0,]$gene
    re = c(length(setdiff(de,he)),length(setdiff(dl,hl)),length(intersect(de,he)),length(intersect(dl,hl)))
    eg_l[[ct]] = setdiff(dl,hl)
    eg_l2[[ct]] = setdiff(de,he)
    }

positive_df = as.data.frame(sapply( eg_l, "[", i = 1:max(sapply( eg_l, length))))
negative_df = as.data.frame(sapply( eg_l2, "[", i = 1:max(sapply( eg_l2, length))))