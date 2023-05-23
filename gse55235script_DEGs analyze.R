Sys.setenv(LANGUAGE="en")
options(stringsAsFactors = F)
rm(list = ls())


library(BiocGenerics)
library(parallel)
library(Biobase)
library(affy)
setwd("E:\\6.22 Hereditas/analyse data and script/DEGs analyze data/GSE55235raw/")
rawdata=ReadAffy()

eset=rma(rawdata)
setwd("E:\\6.22 Hereditas/analyse data and script/DEGs analyze data/raw1 script and data/")
write.exprs(eset,file="NormalizedData.txt")

library(data.table)
eset1=fread("NormalizedData.txt",header = T,data.table = F,sep = "\t")
class(eset1)

BiocManager::install("GEOquery")
library(GEOquery)
setwd("E:\\6.22 Hereditas/analyse data and script/DEGs analyze data/")
rdata.RA1=getGEO(filename = "GSE55235_series_matrix.txt.gz")

ra1=rdata.RA1

pdata1=pData(ra1)
colnames(pdata1)
gpl1=ra1@featureData@data
colnames(gpl1)

gpl1.1=gpl1[,c(11)]
s1=strsplit(gpl1.1,split = " /// ",fixed = T)
s1.1=sapply(s1,function(x){x[1]})
gpl1.2=gpl1[,c(1,11)]
gpl1.2$`Gene Symbol`=s1.1
exp.pl1=merge(gpl1.2,eset1,by=1)
library(dplyr)
exp.pl1.1=distinct(exp.pl1,`Gene Symbol`,.keep_all = T)
exp.pl1.2=na.omit(exp.pl1.1)
row.names(exp.pl1.2)=exp.pl1.2$`Gene Symbol`
exp.pl1.3=exp.pl1.2[,-c(1,2)]

a=c(1:10)
a1=a+20
exp.pl1.4=exp.pl1.3[,c(a,a1)]
group_list1=c(rep("normal",10),rep("RA",10))
group_list1=as.factor(group_list1)
group_list1=relevel(group_list1,ref = "normal")

colnames(exp.pl1.4)=group_list1

library(limma)
exp.pl1.5=normalizeBetweenArrays(exp.pl1.4)
boxplot(exp.pl1.5,outline=FALSE,notch=T,col=group_list1,las=2)

library(RColorBrewer)
colors=brewer.pal(10,"Set3")
boxplot(exp.pl1.5,outline=FALSE,notch=T,col=colors,las=2,ylim=c(4,12))

desgin1=model.matrix(~group_list1)
colnames(desgin1)=levels(group_list1)
row.names(desgin1)=colnames(exp.pl1.5)
fit1=lmFit(exp.pl1.5,desgin1)
fit1.1=eBayes(fit1)
allDiff=topTable(fit1.1,coef = 2,adjust.method = "fdr",number = Inf)

write.csv(allDiff,file = "GSE55235allDiff.csv")

library(ggplot2)
install.packages("ggrepel")
library(ggrepel)
data1=allDiff
data1$significance="stable"

data1$significance[data1$logFC>=1&data1$adj.P.Val<0.05]="up"
data1$significance[data1$logFC<=-1&data1$adj.P.Val<0.05]="down"

p1=ggplot(data1,aes(logFC,-1*log10(adj.P.Val)))+xlim(-6,8)+
  ylim(0,16)+
  geom_point(aes(color=significance),size=5)+
  theme_classic()+
  scale_color_manual(values = c("#00FA9A","#00BFFF","#FF4500"))+
  geom_hline(yintercept = 1.3,linetype=4,size=0.3)+
  geom_vline(xintercept = c(-1,1),linetype=4,size=0.3)+
  theme(title = element_text(size = 22),text = element_text(size = 22))+
  labs(x="log2(foldchange)",y="-log10(P-Value)")

Sys.setenv(LANGUAGE="en")
options(stringsAsFactors = F)
rm(list = ls())

setwd("D:\\bioarticle/6.22/")
load("raw1.RData")

###芯片数据可以不用平均表达量

selcet.log2FC=abs(data1$logFC)>1
table(selcet.log2FC)
select.qval=data1$adj.P.Val<0.05
table(select.qval)

select.vec=(selcet.log2FC&select.qval)
table(select.vec)
dd1=data1[select.vec,]
up1=subset(dd1,dd1$significance=="up")
down1=subset(dd1,dd1$significance=="down")

up1.1=up1[order(up1$adj.P.Val),][1:5,]
down1.1=down1[order(down1$adj.P.Val),][1:5,]

up.down.10=rbind(up1.1,down1.1)

##保存火山图
setwd("D:\\bioarticle/6.22/figure/DEGs/raw1/")
pdf("GSE55235_firemountain.pdf",width = 10,height = 10)
p1+geom_label_repel(data = up.down.10,aes(label=row.names(up.down.10)),
                    box.padding = 1,size=7)
dev.off()

?geom_label_repel

DEGs1=rbind(up1,down1)
write.csv(DEGs1,file = "DEGs1.csv")

annotation_col1=data.frame(Database=c(rep("GEO",20)),
                           CellType=c(rep("normal",10),rep("RA",10)))
exp1.5=exprs(ra1)
exp1.6=exp1.5[,c(a,a1)]

row.names(annotation_col1)=colnames(exp1.6)
class(exp.pl1.5)
exp.pl1.6=as.data.frame(exp.pl1.5)
exprSet.map1=exp.pl1.6[row.names(DEGs1),]
colnames(exprSet.map1)=colnames(exp1.6)

##以矩阵行基因为标准进行scale正态分布
exprSet.map1=t(scale(t(exprSet.map1)))
##为防止分组被聚类隔开，遂用自定义聚类分支顺序，按照分组人为指定热图聚类顺序排序
exprSet.map_t=as.data.frame(t(exprSet.map1))
col_dist=dist(exprSet.map_t)
hclust_1=hclust(col_dist)

dend=reorder(as.dendrogram(hclust_1),
             wts = order(match(rownames(exprSet.map_t),
                               rownames(exprSet.map_t))))
col_cluster=as.hclust(dend)
#自定义分组标注颜色
ann_colors = list(
  CellType = c(normal = "#87CEFA", RA = "#FFB6C1"))

library(pheatmap)
heatmap=pheatmap(exprSet.map1,
         cluster_rows = T,
         cluster_cols = col_cluster,
         annotation_col = annotation_col1,
         annotation_colors = ann_colors,
         show_colnames = F,
         show_rownames = F,
         scale = "row",
         color = colorRampPalette(c("blue","white","red"))(1000))
#pheatmap热图需要自定义函数保存pdf
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
setwd("D:\\bioarticle/6.22/figure/DEGs/raw1/")

save_pheatmap_pdf(heatmap,"GSE55235_pheatmap.pdf",width = 10,height = 10)


###保存矩阵数据、差异分析数据和临床表型数据，用于GSEA和GSVA分析
setwd("D:\\bioarticle/6.22/")
save(data1,exp.pl1.6,exp1.6,file = "GSE55235_GSEA_GSVA_need.Rdata")

##保存矩阵数据用于RBP与DERATGs的相关性分析
setwd("D:\\bioarticle/6.22/RBPDB pertinence analyse data/GSE55235/")
save(exp.pl1.6,exp1.6,file = "RBPDB_pertinence_need.Rdata")

##保存矩阵数据用于CIBERSORT免疫浸润分析
setwd("D:\\bioarticle/6.22/Analysis of immune infiltration in rheumatoid arthritis/CIBERSORT analyse/GSE55235/  ")
save(exp.pl1.6,exp1.6,file = "CIBERSORT_need.Rdata")

##保存矩阵数据用于ssGSEA免疫浸润分析
setwd("D:\\bioarticle/6.22/Analysis of immune infiltration/ssGSEA/GSE55235/")
save(exp.pl1.6,exp1.6,file = "ssGSEA_need.rdata")

#引用R包参考文献
install.packages("pacman")
library(pacman)
p_cite("pheatmap")

