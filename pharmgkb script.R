Sys.setenv(LANGUAGE="en")
options(stringsAsFactors = F)
rm(list = ls())

setwd("D:\\bioarticle/6.22/")
library(data.table)
library(tidyverse)
pharmgkb=unzip("relationships.zip")
setwd("D:\\bioarticle/6.22/PharmGKB gene and drug relationship/")
gene.drug=fread("relationships.tsv",
                header = T,
                data.table = F)

setwd("D:\\bioarticle/6.22/")
?fread
DERATGs=fread("DERATGs.csv",header = T,data.table = F,nrows = Inf)
DERATGs1=data.frame(DERATGs[,-1])
colnames(DERATGs1)="genes"


rm(gene.drug1,gene.drug2)

###用left_join函数筛选gene.drug中和DERATGs匹配的基因信息
?left_join
DERATGs1$score=1
colnames(gene.drug)
colnames(DERATGs1)=c("Entity1_name","scorce")
gene.drug1=left_join(gene.drug,DERATGs1,
                     by="Entity1_name")
gene.drug2=gene.drug1[,c(2,3,5,6,12)]
gene.drug3=gene.drug2%>%
  na.omit()%>%
  subset(Entity2_type=="Chemical")
getwd()
write.csv(gene.drug3,file = "gene_drug.tsv")###保存网络文件



g1=gene.drug3[,1:2]
d1=gene.drug3[,3:4]
colnames(d1)=colnames(g1)
natrue=rbind(g1,d1)
colnames(natrue)=c("name","type")
getwd()
write.csv(natrue,"name_type.tsv")###保存性质文件


##查看基因名称及数目
library(tidyverse)
nam_num_gene=gene.drug3[,-c(2:4)]%>%
  distinct(`Entity1_name`)
setwd("D:\\bioarticle/6.22/analyse data and script/PharmGKB gene and drug data/")
write.csv(nam_num_gene,"nam_num_gene.csv")
##查看药物名称及数目
nam_num_drug=gene.drug3[,-c(1,2,4)]%>%
  distinct(`Entity2_name`)
write.csv(nam_num_drug,"nam_num_drug.csv")
