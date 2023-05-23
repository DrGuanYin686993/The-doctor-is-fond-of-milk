Sys.setenv(LANGUAGE="en")
options(stringsAsFactors = F)
rm(list = ls())

setwd("E:\\6.22 Hereditas/analyse data and script/DrugCentral data/drug.target.interaction.tsv/")
library(data.table)
library(tidyverse)
gene.drug=fread("drug.target.interaction.tsv",
                header = T,
                data.table = F)



DERATGs=read.csv("DERATGs.csv")
DERATGs1=data.frame(DERATGs[,-1])
colnames(DERATGs1)="genes"


###用left_join函数筛选gene.drug中和DERATGs匹配的基因信息
?left_join
DERATGs1$score=1
colnames(gene.drug)
colnames(DERATGs1)=c("GENE","scorce")
gene.drug1=left_join(gene.drug,DERATGs1,
                     by="GENE")
gene.drug2=gene.drug1[,c(1,6,21)]
gene.drug3=gene.drug2%>%
  na.omit()

gene.drug4=gene.drug3%>%
  mutate("feature1"=c("drug"),
         .keep = "all",
         .before = 2)

gene.drug5=gene.drug4%>%
  mutate("feature2"=c("gene"),
         .keep = "all",
         .after = 3)

gene.drug6=gene.drug5[,1:4]

setwd("E:\\6.22 Hereditas/analyse data and script/DrugCentral data/drug.target.interaction.tsv/")
fwrite(gene.drug6,file = "drugcentral_gene and drug.txt")###保存网络文件
