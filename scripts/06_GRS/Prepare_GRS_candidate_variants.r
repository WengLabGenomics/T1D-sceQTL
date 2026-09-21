library(data.table)
geno = data.frame(fread('./input/disease_all_geno/SNP.txt'))
rownames(geno) <- geno$V1
geno = geno[2:ncol(geno)]
remain_snp =read.csv('./res/T1D_risk_eQTL_all.csv')
s_l = c()
for (i in unique(remain_snp$gene)){
    temp = remain_snp[remain_snp$gene == i,]
    temp = temp[order(temp$FDR),]
    s = temp$snps[1]
    s_l = append(s_l,s)
    }
geno = geno[unique(s_l),]

hla = data.frame(fread('./input/HLA_classic/SNP.txt'))
rownames(hla) <- hla$V1
hla = hla[2:ncol(hla)]
remain_allele = read.csv('./res/HLA_condition_res.csv')
hla = hla[remain_allele$Selected_allele,]

GRS_candidate_matrix = rbind(geno,hla)
write.csv(GRS_candidate_matrix,'./grs_input/us_only/Discovery_risk_matrix-us-LD.csv')
