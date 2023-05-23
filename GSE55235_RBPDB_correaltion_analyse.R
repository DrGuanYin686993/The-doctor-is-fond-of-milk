Sys.setenv("LANGUAGE"="en")
options(stringsAsFactors = F)
rm(list = ls())


###读取已经处理过的GEO芯片数据矩阵、筛取的靶基因和RBP数据文件
setwd("D:\\bioarticle/6.22/RBPDB pertinence analyse data/GSE55235/")
load("RBPDB_pertinence_need.Rdata")
colnames(exp.pl1.6)=colnames(exp1.6)

setwd("D:\\bioarticle/6.22/")
library(data.table)
DERATGs=fread("DERATGs1.csv",header = F,data.table = F)

setwd("D:\\bioarticle/6.22/RBPDB pertinence analyse data/RBPDB_csv/")
protein=fread("RBPDB_v1.3.1_proteins_human_2012-11-21.csv",
              header = F,data.table = F)

##处理得到靶基因和RBP各自矩阵
gene1=as.data.frame(t(exp.pl1.6[DERATGs$V1,]))

library(tidyverse)
protein1=protein$V5%>%
  as.data.frame()
colnames(protein1)="Gene Symbol"
protein1$score=1
exp2=exp.pl1.6
exp2$"Gene Symbol"=row.names(exp2)
exp3=left_join(exp2,protein1,by="Gene Symbol")
RBP=exp3%>%
  na.omit()
row.names(RBP)=RBP$`Gene Symbol`
RBP=RBP[,-c(21:22)]
RBP=as.data.frame(t(RBP))
setwd("D:\\bioarticle/6.22/RBPDB pertinence analyse data/GSE55235/")

###RBP与DERATGs相关性分析
library(tidyverse)
cor_result1=list()
#for循环
for (i in 1:ncol(RBP)) {
  exp_cor=cbind(RBP[,i],gene1)
  colnames(exp_cor)[1]=colnames(RBP)[i]
  y=as.numeric(exp_cor[,1])
  colna=colnames(exp_cor)
  cor_data_df=data.frame(colna)
  for (m in 1:length(colna)) {
    test=cor.test(as.numeric(exp_cor[,m]),y,type="speraman")
    cor_data_df[m,2]=test$estimate
    cor_data_df[m,3]=test$p.value
  }
  names(cor_data_df)=c("symbol","correlation","pvalue")
  cor_result1[[i]]=cor_data_df[,-3]
  print(i)
}

cor_result1[[1]]


gene_RBP=data.frame(rep("1",54))
row.names(gene_RBP)=DERATGs$V1
for (i in 1:length(cor_result1)) {
  a=as.data.frame(cor_result1[[i]][-1,])
  row.names(a)=a$symbol
  a=as.data.frame(a[,-1])
  colnames(a)=cor_result1[[i]]$symbol[1]
  row.names(a)=cor_result1[[i]]$symbol[-1]
  gene_RBP[,i+1]=a
  row.names(gene_RBP)=cor_result1[[i]]$symbol[-1]
}

gene_RBP=gene_RBP[,-1]
gene_RBP1=as.data.frame(t(gene_RBP))

###画热图
library(pheatmap)
RBP_DERATGs_cor=pheatmap(gene_RBP1,
             show_colnames = T,
             show_rownames = T,
             scale = "row",
             cluster_rowsu = T,
             cluster_cols = T,
             colorRampPalette(c("#000080","#AFEEEE",
                                "#FAF0E6","#FFA07A",
                                "#FF4500"))(100),
             fontsize=20)

#pheatmap热图需要自定义函数保存pdf
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
setwd("D:\\bioarticle/6.22/figure/RBP correlation analyse/")

save_pheatmap_pdf(RBP_DERATGs_cor,"GSE55235_RBP_pheatmap.pdf",width = 20,height = 80)


setwd("D:\\bioarticle/6.22/RBPDB pertinence analyse data/GSE55235/")
