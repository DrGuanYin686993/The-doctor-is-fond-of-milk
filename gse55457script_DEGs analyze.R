Sys.setenv("LANGUAGE"="en")
options(stringsAsFactors = F)
rm(list = ls())

setwd("E:\\6.22 Hereditas/analyse data and script/DEGs analyze data/raw2 script and data/")

library(BiocGenerics)
library(parallel)
library(Biobase)
library(affy)
rawdata=ReadAffy()

eset=rma(rawdata)
write.exprs(eset,file="raw2 NormalizedData.txt")
library(data.table)
exp.1=fread("raw2 NormalizedData",header = T,data.table = F)

?fread

setwd("E:\\6.22 Hereditas/analyse data and script/DEGs analyze data/")
library(GEOquery)
rdata.RA2=getGEO(filename = "GSE55457_series_matrix.txt.gz")

ra2=rdata.RA2
exp1=exprs(ra2)
load("raw2.RData")
pdata=pData(ra2)

colnames(pdata)

gpl=ra2@featureData@data
colnames(gpl)
gpl1=gpl[,c(1,11)]

a=strsplit(gpl1$`Gene Symbol`,split =" /// ",fixed = T )
a1=sapply(a,function(x){x[1]})
gpl1$`Gene Symbol`=a1

exp.pl.1=merge(gpl1,exp.1,by=1)
library(dplyr)
exp.pl.2=distinct(exp.pl.1,`Gene Symbol`,.keep_all = T)
exp.pl.3=na.omit(exp.pl.2)
row.names(exp.pl.3)=exp.pl.3$`Gene Symbol`
exp.pl.4=exp.pl.3[,-c(1,2)]
b=c(1:10)
b1=c(11:23)
colnames(exp.pl.4)=colnames(exp1)
exp.pl.5=exp.pl.4[,c(b,b1)]

group_list=c(rep("normal",10),rep("RA",13))
group_list1=as.factor(group_list)
group_list2=relevel(group_list1,ref = "normal")

library(limma)
exp.pl.6=normalizeBetweenArrays(exp.pl.5)
library(RColorBrewer)
colors=brewer.pal(10,"Set3")

boxplot(exp.pl.6,outline=FALSE,notch=T,col=colors,
        las=2)

design=model.matrix(~group_list2)
colnames(design)=levels(group_list2)
row.names(design)=colnames(exp.pl.6)
fit=lmFit(exp.pl.6,design)
fit=eBayes(fit)
allDiff2=topTable(fit,coef = 2,adjust.method = "fdr",number = Inf)

write.csv(allDiff2,file = "GSE55457allDiff.csv")

library(ggplot2)
library(ggrepel)
data=allDiff2

data$significance="stable"
data$significance[data$logFC>=1&data$adj.P.Val<0.05]="up"
data$significance[data$logFC<=-1&data$adj.P.Val<0.05]="down"

p=ggplot(data,aes(logFC,-1*log10(adj.P.Val)))+
  xlim(-4,6)+
  ylim(0,5)+
  geom_point(aes(color=significance),size=5)+
  theme_classic()+
  scale_color_manual(values = c("#00FA9A","#00BFFF","#FF4500"))+
  geom_hline(yintercept = 1.3,linetype=4,size=0.3)+
  geom_vline(xintercept = c(-1,1),linetype=4,size=0.3)+
  theme(title = element_text(size = 22),text = element_text(size = 22))+
  labs(x="log2(foldchange)",y="-log10(P.Value)")

select.log2FC=abs(data$logFC)>1
select.qval=data$adj.P.Val<0.05
select.vec=(select.log2FC&select.qval)
table(select.vec)

dd1=data[select.vec,]
write.csv(dd1,file = "DEGs2.csv")
up=subset(dd1,significance=="up")
donw=subset(dd1,significance=="down")
up1=up[order(dd1$adj.P.Val),][1:5,]
down1=donw[order(dd1$adj.P.Val),][1:5,]
up.down=rbind(up1,down1)

##保存火山图
setwd("D:\\bioarticle/6.22/figure/DEGs/raw2/")
pdf("GSE55457_firemountain.pdf",width = 10,height = 10)
p+geom_label_repel(data = up.down,aes(label=row.names(up.down)),
                   box.padding = 1,size=7)
dev.off()



annotation_col1=data.frame(Database=c(rep("GEO",23)),
                           CellType=c(rep("normal",10),
                                      rep("RA",13)))
row.names(annotation_col1)=colnames(exp.pl.6)
exprSet.map=exp.pl.6[row.names(dd1),]
library(pheatmap)
##以矩阵行基因为标准进行scale正态分布
exprSet.map=t(scale(t(exprSet.map)))
##由于分组被聚类隔开，遂用自定义聚类分支顺序，按照分组人为指定热图聚类顺序排序
exprSet.map_t=as.data.frame(t(exprSet.map))
col_dist=dist(exprSet.map_t)
hclust_1=hclust(col_dist)

dend=reorder(as.dendrogram(hclust_1),
             wts = order(match(rownames(exprSet.map_t),
                               rownames(exprSet.map_t))))
col_cluster=as.hclust(dend)
#自定义分组标注颜色
ann_colors = list(
  CellType = c(normal = "#87CEFA", RA = "#FFB6C1"))

heatmap=pheatmap(exprSet.map,
         cluster_rows = T,
         cluster_cols = col_cluster,
         annotation_col=annotation_col1,
         annotation_colors = ann_colors,
         show_colnames=F,
         show_rownames = F,
         color=colorRampPalette(c("blue","white","red"))(100))
?pheatmap

#pheatmap热图需要自定义函数保存pdf
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
setwd("D:\\bioarticle/6.22/figure/DEGs/raw2/")

save_pheatmap_pdf(heatmap,"GSE55457_pheatmap.pdf",width = 10,height = 10)


###保存矩阵数据、差异分析数据和临床表型数据，用于GSEA和GSVA分析
setwd("D:\\bioarticle/6.22/")
save(allDiff2,exp.pl.6,pdata,file = "GSE55457_GSEA_GSVA_need.rdata")
getwd()

##保存exp.pl.6矩阵数据用于RBP与DERATGs的相关性分析
setwd("D:\\bioarticle/6.22/RBPDB pertinence analyse data/GSE55457/")
save(exp.pl.6,file = "GSE55457_RBP_pertinence_need.rdata")

##保存exp.pl.6矩阵数据、pdata临床表型数据，用于CIBERSORT免疫浸润分析
setwd("D:\\bioarticle/6.22/Analysis of immune infiltration/CIBERSORT analyse/GSE55457/")
save(exp.pl.6,pdata,file = "CIBERSORT_need.rdata")

##保存矩阵数据用于ssGSEA免疫浸润分析
setwd("D:\\bioarticle/6.22/Analysis of immune infiltration/ssGSEA/GSE55457/")
save(exp.pl.6,pdata,file = "ssGSEA_need.rdata")
