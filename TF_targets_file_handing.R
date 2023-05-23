Sys.setenv(LANGUAGE="en")
options(stringsAsFactors = F)
rm(list = ls())
setwd("D:\\bioarticle/6.22/TRRUST TF-targets data/")

library(data.table)
library(tidyverse)
TF.targets=fread("c4742ec51a.key_regulators.tsv",header = T,data.table = F)
tf1=TF.targets[,c(1,6)]

judge_na_num=function(x)sum(is.na(x))!=39
largest_gene_num=1000
judge_na_num(tf1$`List of overlapped genes`)
judge_na_num(tf1$`Key TF`)
colnames(tf2)
tf2=tf1%>%
  separate(col = `List of overlapped genes`,
           into = paste0("x",1:largest_gene_num),
           sep = ",")%>%
  select(where(judge_na_num))%>%
  pivot_longer(num_range("x",1:largest_gene_num),
               values_to = "node2")%>%
  na.omit()%>%
  select(!name)%>%
  rename("node1"="Key TF")%>%
  group_by(node1,node2)%>%
  summarise(interact=n())

getwd()  
write.csv(tf2,"tf_tarfets_cyto.tsv")

#查看DERATGs基因名称和数目
str(tf2)
node2_gene=tf2[,-1]
node2_gene_num=distinct(node2_gene,`node2`)
setwd("D:\\bioarticle/6.22/analyse data and script/TRRUST TF-targets data/")
write.csv(node2_gene_num,"DERATGs_node2_genes.csv")
