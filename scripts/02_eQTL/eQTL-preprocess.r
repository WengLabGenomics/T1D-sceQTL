library(dplyr)
library(data.table)
library(stringr)

ref_clean = read.csv('./input/gene_ref_clean.csv',row.names = 1)
write.table(ref_clean,'./input/disease_all_geno/geneloc.txt',quote = FALSE)

geno = data.frame(fread('./input/snp_qc_all.csv',header = TRUE))
rownames(geno) =geno$V1
geno = geno[,2:ncol(geno)]
colnames(geno) = sub(".", "", colnames(geno))
write.table(geno,'./input/ddisease_all_geno/SNP.txt',quote = FALSE)

first_part <- sapply(strsplit(rownames(geno), ":"), function(x) x[1])
sec_part <- sapply(strsplit(rownames(geno), ":"), function(x) x[2])
snploc = data.frame(row.names = rownames(geno),'chr' = first_part,'pos' = sec_part)
snploc$chr = paste0('chr',snploc$chr)
snploc = data.frame(row.names = rownames(geno),'chr' = first_part,'pos' = sec_part)
snploc$pos = as.numeric(snploc$pos)
snploc$chr = paste0('chr',snploc$chr)
                   
write.table(snploc,'./input/disease_all_geno/snpsloc.txt',quote = FALSE)
                   
cli  = read.csv('./input//all_cli_v3.csv')
rownames(cli) <- cli$RNA_id
cli$sex = as.numeric(cli$sex)
cli$age = as.numeric(cli$age)
cli$batch = as.numeric(cli$batch)
cli[cli$Disease == 'T1D','Disease'] = '1'
cli[cli$Disease == 'HC','Disease'] = '0'
cli$Disease = as.numeric(cli$Disease)                   

#geno pca
pca = read.table('./input/snp_qc_0.8.eigenvec')
rownames(pca) <- pca$V2
pca = pca[,3:ncol(pca)]
pca = data.frame(t(pca))
cli2 = cli[!is.na(cli$geno_map_id),]
rownames(cli2) = cli2$geno_map_id
colnames(pca) = as.character(cli2[colnames(pca),'RNA_id'])
rownames(pca) = paste0(rep('genoPC',nrow(pca)),seq(1,nrow(pca)))
cli3 = cli[c('sex','age','batch','Disease')]
cli3 = data.frame(t(cli3))                   
colnames(pca) <- paste0('X',colnames(pca))
pca = pca[colnames(cli3)]
cli4 = rbind(cli3,pca[c('genoPC1','genoPC2','genoPC3'),])
cli4 = cli4[c('sex','age','batch','genoPC1','genoPC2','genoPC3','Disease'),]                   
colnames(cli4) = sub(".", "", colnames(cli4))                   
write.table(cli4,'./input/disease_all_geno/Covariates.txt',quote = FALSE)                   

                   
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
data = data.frame(t(data))
colnames(data) = sub(".", "", colnames(data)) 
write.table(data,paste0('./input/disease_all_geno/',ct,'_GE.txt'),quote = FALSE)
    }  
