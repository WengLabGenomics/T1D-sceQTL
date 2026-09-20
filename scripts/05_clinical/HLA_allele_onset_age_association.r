library(stringr)
library(data.table)
library(mice)
library(RNOmni)
library(vegan)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(Hmisc)

merge_df <- read.csv('./HLA_matrix/Discovery_cohor_HLA_dosage.csv',row.names = 1)
merge_df[merge_df$Disease == 'T1D','disease_status'] = 1
merge_df[merge_df$Disease == 'HC','disease_status'] = 0
merge_df$disease_status = as.numeric(merge_df$disease_status)
test = merge_df[grepl('HLA',colnames(merge_df))]
test$sex = merge_df$sex
test$onset_age = merge_df$onset_age
test$province = merge_df2$province
test$onset_age = merge_df$onset_age

name = c()
p <- c()
coe <- c()
OR <- c()
CI_L = c()
CI_H = c()
SE = c()
for( i in seq(1:(length(test)-length(c('sex','onset_age','province'))))){
    temp = test[c(colnames(test)[i],'sex','onset_age','province')]
    colnames(temp) <- c('loci','sex','onset_age','province')
    temp = temp[!is.na(temp['onset_age']),]
    mean_value <- mean(temp[,'onset_age'], na.rm = TRUE)
    sd_value <- sd(temp[,'onset_age'], na.rm = TRUE)
    lower_bound <- mean_value - 5 * sd_value
    upper_bound <- mean_value + 5 * sd_value
    temp = temp[temp['onset_age'] > lower_bound &  temp['onset_age'] <upper_bound, ]
    temp['onset_age'] =  RankNorm(temp[,'onset_age'])
    model <- lm(onset_age ~loci+sex+province, data = temp)
    null_model <- lm(onset_age ~ sex+province, data = temp)
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
rownames(HLA_res) <-HLA_res$name
write.csv(HLA_res,'./HLA_res/T1D_onset_age_HLA_association_res.csv')

can  <- read.csv('./HLA_res/disease_sig_hla_loci.csv',row.names = 1)
can = can[can$p < 5e-8,]
rownames(can) <- can$name

os = HLA_res[rownames(can),]

data <- data.frame(
    row.names = rownames(can),
  onset = os$coe,
  T1D = can$coe,
  group = can$sig2,
    se_x = os$SE,
    se_y = can$SE,
  error_x_L = os$CI_L,
    error_x_H = os$CI_H,
  error_y_L = can$CI_L,
    error_y_H = can$CI_L
)
library(stringr)
x1 = str_split_fixed(rownames(data),'_',3)
x2 = str_split_fixed(x1[,3],'\\.',2)
re = sapply(x2[,1], function(x) {
 if (grepl("^\\d{2,}$", x)) {
    return(x)  
  } else {
    return(paste0("0", x)) 
  }
})
names(re) <- NULL
data$name = paste0(x1[,2],'*',re,':',x2[,2])

data[rownames(data) %in% setdiff(rownames(can),HLA_res[!is.na(HLA_res$p) & HLA_res$p < 0.05 & HLA_res$name %in% rownames(can),]$name),'group'] = 'non-sig'
data[data$group %in% c('susceptive','protective'),'label'] = data[data$group %in% c('susceptive','protective'),'name']


p <- ggplot(data, aes(x = onset, y = T1D, color = group)) +
  geom_point(size = 3) +  
  geom_errorbar(aes(xmin =  onset-(1.96*se_x), xmax = onset+(1.96*se_x)), width = 0.02) + 
  geom_errorbar(aes(ymin =  T1D-(1.96*se_y), ymax = T1D+(1.96*se_y)), height = 0.02) + 
  theme_bw() +theme(panel.grid = element_blank(),aspect.ratio = 1)+
  labs(x = '', y = '') +
  theme(legend.position = "top")+
geom_vline(xintercept = 0, linetype = "dashed", color = "gray")+
geom_hline(yintercept = 0, linetype = "dashed", color = "gray")+
scale_color_manual(values = c("susceptive" = '#D6604D', "protective" = '#4393C3','non-sig' ='lightgrey'))+
geom_text_repel(aes(label = label),color = 'black')
p


merge_df <- read.csv('./HLA_matrix/Discovery_cohor_HLA_dosage.csv',row.names = 1)
T1D = merge_df[merge_df$Disease == 'T1D',]
can  <- read.csv('./HLA_res/disease_sig_hla_loci.csv',row.names = 1)
can = can[can$p < 5e-8,]
can = can[can$sig2 == 'susceptive',]
df = data.frame(rowSums(T1D[can$name]))
colnames(df) = 'Counts'
df[c('onset_age','age','sex','HbA1c','IA2','ZnT8','GAD','Duration')] = T1D[rownames(df),c('onset_age','age','sex','HbA1c','IA2','ZnT8','GAD','Duration')]
res = rcorr(as.matrix(df),type = 'spearman')

library(corrplot)
corrplot(t(r), p.mat = t(p), method = 'color', diag = FALSE,
         sig.level = 0.01,pch.cex = 0.9, insig = 'label_sig', pch.col = 'black',tl.col = "black",col = rev(COL2('RdBu', 200)))
