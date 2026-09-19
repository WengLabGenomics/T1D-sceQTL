library(stringr)
library(data.table)
library(mice)
library(RNOmni)
library(vegan)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)

df <- read.csv('./HLA_matrix/Discovery_cohor_HLA_matrix.csv',row.names = 1)
T1D = df[df$Disease  == 'T1D',]

shannon_df = data.frame(row.names  = c('HLA','shannon_index'))
for (i in c('A','B','C','DRB1','DQA1','DQB1','DPA1','DPB1')){#c('A','B','C','DRB1','DQA1','DQB1','DPA1','DPB1')
i1 = paste0(i,'.Allele1')
i2 = paste0(i,'.Allele2')
test = T1D[c(i1,'disease_status','sex')]
colnames(test) <- c('HLA','disease_status','sex')
vec = test$HLA
vec <- ifelse(substr(vec, 1, 1) == "0", substring(vec, 2), vec)
test$HLA = vec
re = c()
for (x in seq(1:length(test$HLA))){
 l = length(strsplit(test$HLA, ':')[[x]])
    re_x = ifelse(l>2,paste0(strsplit(test$HLA, ':')[[x]][1],':',strsplit(test$HLA, ':')[[x]][2]),
                 test$HLA[x])
    re = append(re,re_x)
}
    test$HLA = re

dummy_data <- model.matrix(~0+ HLA, data = test)
dummy_data = data.frame(dummy_data)
test = T1D[c(i2,'disease_status','sex')]
colnames(test) <- c('HLA','disease_status','sex')
vec = test$HLA
vec <- ifelse(substr(vec, 1, 1) == "0", substring(vec, 2), vec)
test$HLA = vec

        re = c()
for (x in seq(1:length(test$HLA))){
 l = length(strsplit(test$HLA, ':')[[x]])
    re_x = ifelse(l>2,paste0(strsplit(test$HLA, ':')[[x]][1],':',strsplit(test$HLA, ':')[[x]][2]),
                 test$HLA[x])
    re = append(re,re_x)
}
    test$HLA = re

dummy_data2 <- model.matrix(~0+ HLA, data = test)
dummy_data2 = data.frame(dummy_data2)
v1 <- intersect(colnames(dummy_data),colnames(dummy_data2))
v2 <- colnames(dummy_data[, colSums(dummy_data)>= 0])
v3 <- colnames(dummy_data2[, colSums(dummy_data2)>=0])
inte = dummy_data[v1]+dummy_data2[v1]
a1 =dummy_data[setdiff(v2,v1)]
a2=dummy_data2[setdiff(v3,v1)]
    if (dim(a1)[2] > 0 & dim(a2)[2] > 0 ){
        test = cbind(inte,c(a1,a2))}
    else if (dim(a1)[2] > 0 & dim(a2)[2]  ==  0 ){
        test = cbind(inte,a1)}
    else if (dim(a1)[2] == 0 & dim(a2)[2]  >  0 ){
        test = cbind(inte,a2)}
    else{
        test = inte}
shannon_index <- diversity(data.frame(colSums(test))[,1], index = "shannon")
shannon_df_temp = data.frame('HLA' = i,'shannon_index' = shannon_index)
    if (i == 'A'){
       shannon_df =  shannon_df_temp
        }else{shannon_df = rbind(shannon_df,shannon_df_temp)}
    }
d_t1d =  shannon_df
d_t1d$disease = 'T1D'

HC = df[df$Disease  == 'HC',]
shannon_df = data.frame(row.names  = c('HLA','shannon_index'))
for (i in c('A','B','C','DRB1','DQA1','DQB1','DPA1','DPB1')){
i1 = paste0(i,'.Allele1')
i2 = paste0(i,'.Allele2')
test = HC[c(i1,'disease_status','sex')]
colnames(test) <- c('HLA','disease_status','sex')
vec = test$HLA
vec <- ifelse(substr(vec, 1, 1) == "0", substring(vec, 2), vec)
test$HLA = vec

re = c()
for (x in seq(1:length(test$HLA))){
 l = length(strsplit(test$HLA, ':')[[x]])
    re_x = ifelse(l>2,paste0(strsplit(test$HLA, ':')[[x]][1],':',strsplit(test$HLA, ':')[[x]][2]),
                 test$HLA[x])
    re = append(re,re_x)
}
    test$HLA = re

dummy_data <- model.matrix(~0+ HLA, data = test)
dummy_data = data.frame(dummy_data)
test = HC[c(i2,'disease_status','sex')]
colnames(test) <- c('HLA','disease_status','sex')
vec = test$HLA
vec <- ifelse(substr(vec, 1, 1) == "0", substring(vec, 2), vec)
test$HLA = vec

        re = c()
for (x in seq(1:length(test$HLA))){
 l = length(strsplit(test$HLA, ':')[[x]])
    re_x = ifelse(l>2,paste0(strsplit(test$HLA, ':')[[x]][1],':',strsplit(test$HLA, ':')[[x]][2]),
                 test$HLA[x])
    re = append(re,re_x)
}
    test$HLA = re

dummy_data2 <- model.matrix(~0+ HLA, data = test)
dummy_data2 = data.frame(dummy_data2)
v1 <- intersect(colnames(dummy_data),colnames(dummy_data2))
v2 <- colnames(dummy_data[, colSums(dummy_data)>= 0])
v3 <- colnames(dummy_data2[, colSums(dummy_data2)>=0])
inte = dummy_data[v1]+dummy_data2[v1]
a1 =dummy_data[setdiff(v2,v1)]
a2=dummy_data2[setdiff(v3,v1)]
    if (dim(a1)[2] > 0 & dim(a2)[2] > 0 ){
        test = cbind(inte,c(a1,a2))}
    else if (dim(a1)[2] > 0 & dim(a2)[2]  ==  0 ){
        test = cbind(inte,a1)}
    else if (dim(a1)[2] == 0 & dim(a2)[2]  >  0 ){
        test = cbind(inte,a2)}
    else{
        test = inte}
shannon_index <- diversity(data.frame(colSums(test))[,1], index = "shannon")
shannon_df_temp = data.frame('HLA' = i,'shannon_index' = shannon_index)
    if (i == 'A'){
       shannon_df =  shannon_df_temp
        }else{shannon_df = rbind(shannon_df,shannon_df_temp)}
    }
d_hc = shannon_df
d_hc$disease = 'HC'

d = rbind(d_t1d,d_hc)
d$HLA = factor(d$HLA,levels = unique(d$HLA))
color_m = c('#FB8072','#FDB462','#FCCDE5','#8DD3C7','#BEBADA','#80B1D3','#CCEBC5','#B3DE69')
p1 <- ggplot(d, aes(x = disease,y = shannon_index)) +
  geom_boxplot() + theme_bw()+
  geom_dotplot(
    aes(fill = HLA), trim = FALSE,
    binaxis='y', stackdir='centerwhole', dotsize = 1.2,
    position = position_dodge(0.2))+
scale_fill_manual(values=color_m)+
labs(title = "", x = "Number of digits of HLA alleles", y = "Shannon's Diversity Index")+
geom_line(aes(group = HLA), color = 'gray', lwd = 0.5)+

theme(aspect.ratio = 2,
        axis.line = element_line(color = "black",size = 0.4),
        axis.text.y = element_text(color="black",size=10),
        axis.text.x = element_text(margin = margin(t = 2),color="black",size=10),
        panel.spacing = unit(0,"lines"),
     panel.grid = element_blank())+ 
  stat_compare_means(method = "t.test",paired = TRUE, 
                     comparisons=list(c("T1D", "HC")))
p1