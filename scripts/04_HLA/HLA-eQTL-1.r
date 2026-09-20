library(dplyr)
library(data.table)
library(stringr)

merge_df <- read.csv('./HLA_matrix/Discovery_cohor_HLA_dosage.csv',row.names = 1)
geno =  merge_df[grepl('HLA',colnames(merge_df))]
write.table(geno,'./input/HLA_classic/SNP.txt',quote = FALSE)


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
write.table(cli4,'./input/HLA_classic/Covariates.txt',quote = FALSE)                 

                   
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
write.table(data,paste0('./input/HLA_classic/',ct,'_GE.txt'),quote = FALSE)
    }        

loci = c('A','B','C','DRB1','DQA1','DQB1','DPB1','DPA1')
pos = c(29944214,31355516,31270415,32584287,32640529,32663284,33081591,33069309)
df = data.frame('loci'=loci,pos = pos)
x= data.frame(str_split_fixed(rownames(geno),'_',3))
colnames(x) <- c('X1','loci','X3')
rownames(x) = rownames(geno)

x = left_join(x,df,by = 'loci')
rownames(x) = rownames(geno)

x$chr = 'chr6'

snploc = x[c('chr','pos')]

write.table(snploc,'./input/HLA_classic/snpsloc.txt',quote = FALSE)

ref_clean = read.csv('./input/gene_ref_clean.csv',row.names = 1)
write.table(ref_clean,'./input/HLA_classic/geneloc.txt',quote = FALSE)


library(MatrixEQTL)

for (ct in c('CD8T', 'NK', 'CD4T', 'CD14Mono', 'B', 'CD16Mono', 'pDC', 'cDC','T_Prolif', 'Plasma')) {
base.dir = './input//HLA_classic'
save.dir = './res/HLA_classic_All'
useModel = modelLINEAR

# Genotype file name
SNP_file_name = paste0(base.dir,"/SNP.txt");
snps_location_file_name = paste0(base.dir, "/snpsloc.txt");

# Gene expression file name
expression_file_name = paste0(base.dir,'/',ct, "_GE.txt");
gene_location_file_name = paste0(base.dir, "/geneloc.txt");

# Covariates file name
# Set to character() for no covariates
covariates_file_name = paste0(base.dir,"/Covariates.txt");

# Output file name
output_file_name_cis = tempfile();
output_file_name_tra = tempfile();

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 1;
pvOutputThreshold_tra = 1e-2;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

# Distance for local gene-SNP pairs
cisDist = 1e6;

## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = " ";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = " ";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);


## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = " ";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
    cvrt$LoadFile(covariates_file_name);
}

## Run the analysis
snpspos = read.table(snps_location_file_name,
                        header = TRUE, stringsAsFactors = FALSE);
snpspos$snpid = rownames(snpspos)
snpspos = snpspos[c('snpid','chr','pos')]
genepos = read.table(gene_location_file_name,
                        header = TRUE, stringsAsFactors = FALSE);
genepos$geneid = rownames(genepos)
genepos = genepos[c('geneid','chr','left','right')]

# 提取样本名（列名）
snp_samples <- snps$columnNames
gene_samples <- gene$columnNames
common_samples <- intersect(snp_samples, gene_samples)

# 筛选 snps 样本
snps$ColumnSubsample(match(common_samples,snps$columnNames))

# 筛选 gene 样本
gene$ColumnSubsample(match(common_samples,gene$columnNames))

cvrt$ColumnSubsample(match(common_samples,cvrt$columnNames))

me = Matrix_eQTL_main(
        snps = snps,
        gene = gene,
        cvrt = cvrt,
        output_file_name  = output_file_name_tra,
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
        noFDRsaveMemory = FALSE);


saveRDS(me,paste0(save.dir,'/',ct, "_eQTL.rds"))

res = me$cis$eqtls
res = res[res$FDR < 0.05,]
write.csv(res,paste0(save.dir,'/',ct, "_cis_fdr005.csv"))
    }