library(stringr)
library(data.table)
library(mice)
library(RNOmni)
library(vegan)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)

merge_df <- read.csv('./HLA_matrix/Discovery_cohor_HLA_dosage.csv',row.names = 1)
merge_df[merge_df$Disease == 'T1D','disease_status'] = 1
merge_df[merge_df$Disease == 'HC','disease_status'] = 0
merge_df$disease_status = as.numeric(merge_df$disease_status)

test = merge_df[grepl('HLA',colnames(merge_df))]
test = test[,colSums(test)>= (nrow(test)*2*.01)]
test$sex = merge_df$sex
test$disease_status = merge_df$disease_status
test$province = merge_df$province

name = c()
p <- c()
coe <- c()
OR <- c()
CI_L = c()
CI_H = c()
SE = c()
for( i in seq(1:(length(test)-length(c('sex','disease_status','province'))))){
#for( i in seq(1:4)){
    temp = test[c(colnames(test)[i],'sex','disease_status','province')]
    colnames(temp) <- c('loci','sex','disease_status','province')
    model <- glm(disease_status ~loci+sex+province, family = binomial(link = "logit"), data = temp)
    null_model <- glm(disease_status ~ sex+province, family = binomial(link = "logit"), data = temp)
    likelihood_ratio_test <- anova(null_model,model, test = "Chisq")
    p_value <- likelihood_ratio_test$Pr[2]
    or <- exp(coef(model))['loci']
    se <- sqrt(diag(vcov(model)))['loci']
    ci_low =exp(confint(model, "loci"))[1]
    ci_high =exp(confint(model, "loci"))[2]
    name = append(name,colnames(test)[i])
    p <- append(p,p_value)
    OR <- append(OR,or)
    coe <- append(coe,coef(model)["loci"])
    CI_L = append(CI_L,ci_low)
    CI_H = append(CI_H,ci_high)
    SE = append(SE,se)
    }
HLA_res = data.frame(name,p,coe,OR,CI_L,CI_H,SE)


library(stringr)
x1 = str_split_fixed(HLA_res$name,'_',3)
x2 = str_split_fixed(x1[,3],'\\.',2)
re = sapply(x2[,1], function(x) {
 if (grepl("^\\d{2,}$", x)) {
    return(x)  
  } else {
    return(paste0("0", x)) 
  }
})
names(re) <- NULL
HLA_res$name2 = paste0(x1[,2],'*',re,':',x2[,2])

HLA_res_all = HLA_res
HLA_res_all$loci = str_split_fixed(HLA_res$name,'_',3)[,2]
HLA_res_all$logP = -log10(HLA_res_all$p)
HLA_res_all[HLA_res_all$coe >0 & HLA_res_all$p < 5e-8,'sig2'] = 'susceptive'
HLA_res_all[HLA_res_all$coe <0 & HLA_res_all$p < 5e-8,'sig2'] = 'protective'
res1 = HLA_res_all
HLA_res_all = res1
HLA_res_all$loci <- factor(HLA_res_all$loci,levels = c('A','B','C','DRB1','DQA1','DQB1','DPA1','DPB1'))
HLA_res_all$sig2 <- factor(HLA_res_all$sig2,levels = c('susceptive','protective'))
HLA_res_all[!is.na(HLA_res_all$sig2),'label'] = HLA_res_all[!is.na(HLA_res_all$sig2),'name2']
HLA_res_all[HLA_res_all$coe < 0,'logP'] = -HLA_res_all[HLA_res_all$coe < 0,'logP']

library(ggplot2)
library(RColorBrewer)
library(ggrepel)
p1 <- ggplot(data = HLA_res_all[HLA_res_all$coe > 0,],aes(x=loci,y=logP,colour = sig2))+
geom_point(size = 3.5,alpha  = 1,position = position_jitter(width = 0.4, height = 0.05))+theme_classic()+
scale_color_manual(values = c("susceptive" = '#D6604D', "protective" = '#4393C3'))+
labs(color = "", x = '',y = '-log10(P-value)')+
theme(panel.grid = element_blank(), axis.ticks.x = element_line(color = "black"),legend.position = "none")+#ylim(-0.5,12.5)+
geom_text_repel(aes(label = label),color = 'black')+
geom_hline(yintercept = c(-log10(5e-8)), linetype = "dashed", color = "black", linewidth = 0.5)
p2 <- ggplot(data = HLA_res_all[HLA_res_all$coe < 0,],aes(x=loci,y=logP,colour = sig2))+
geom_point(size = 3.5,alpha  = 1,position = position_jitter(width = 0.4, height = 0.05))+theme_classic()+
scale_color_manual(values = c("susceptive" = '#D6604D', "protective" = '#4393C3'))+
labs(color = "", x = '',y = '-log10(P-value)')+
theme(panel.grid = element_blank(), axis.ticks.x = element_line(color = "black"),legend.position = "none")+#ylim(-12.5,0.5)+
geom_text_repel(aes(label = label),color = 'black')+
geom_hline(yintercept = log10(5e-8), linetype = "dashed", color = "black", linewidth = 0.5)
library(patchwork)
p1/p2

write.csv(HLA_res_all,'./HLA_res/disease_sig_hla_loci.csv')