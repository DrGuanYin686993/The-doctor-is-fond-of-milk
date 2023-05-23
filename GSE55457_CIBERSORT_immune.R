Sys.setenv(LANGUAGE="en")
options(stringsAsFactors = F)
rm(list = ls())


###读取预先处理过的GSE55457数据集（raw原始数据）的
###所需矩阵数据，以及临床表型数据
setwd("D:\\bioarticle/6.22/Analysis of immune infiltration/CIBERSORT analyse/GSE55457/")
load("CIBERSORT_need.rdata")

###矩阵处理
###处理要求：表达矩阵格式第一列是基因名，第一行是样品名，
###不能有重复基因名，第一列列名不能有空白。
###矩阵中不能存在空白或NA值，不要对表达量取log2
##之前处理的是原始数据raw，并且已经用rma（）标准化过了
##注：rma使用pm信号，会对exp数据进行log2处理，
##所以此处需要对log2处理过的矩阵进行2^返回处理
library(tidyverse)
exp1=2^exp.pl.6%>%
  as.data.frame()

setwd("D:\\bioarticle/6.22/Analysis of immune infiltration/CIBERSORT analyse/GSE55457/three files/")
??fwrite
library(data.table)
fwrite(exp1,"DATA.txt",sep = "\t",
       col.names = T,row.names = T)

DATA=fread("DATA.txt",header = T,data.table = F,
           sep = "\t")###打开保存的exp1矩阵看是否符合矩阵处理要求，发现不符合则重新处理
DATA=rename(DATA,"Gene Symbol"="V1")
fwrite(DATA,file = "DATA.txt",sep = "\t",
       row.names = F,col.names = T)###处理完原exp1数据再次保存为原文件名文件

###LM22.txt文件已在GSE55235中脚本获取，直接文件夹复制粘贴即可，遂此步具体见GSE55235分析脚本
#打开从nature文章下载的参考数据集（该参考数据集，提供22种免疫细胞亚型的基因表达特征集：LM22。
#LM22.txt获取方法：
#在Cibersort论文中下载Supplementry table 1
#https://www.nature.com/articles/nmeth.3337#MOESM207
#将LM22处理成满足上述矩阵处理要求的矩阵，并保存为原文件名文件

###Cibersort.R已在GSE55235获取，直接复制粘贴即可
#在R中新建R Script，复制以下网址中代码，保存为“Cibersort.R”
#https://rdrr.io/github/singha53/amritr/src/R/supportFunc_cibersort.R

###将”DATA.txt",“LM22.txt","Cibersort.R"保存在同一文件夹，运行Cibersort的代码如下：
setwd("D:\\bioarticle/6.22/Analysis of immune infiltration/CIBERSORT analyse/GSE55457/three files/")
source("Cibersort.R")
result=CIBERSORT("LM22.txt","DATA.txt",perm = 1000,QN=T)
#perm置换次数=1000，QN分位数归一化=TRUE                
#QN如果是芯片设置为T，如果是测序就设置为F  

###注意：Cibersort结果的默认文件名为CIBERSORT-Results.txt，
###在同一文件夹下进行第二次运算会覆盖第一次得到的文件，
###建议在每一次运算之后对文件重命名。

####ggheatmap免疫浸润热图展示
###预先处理得到所需数据格式（矩阵）
exp_result1=as.matrix(t(result[,1:22]))

##以矩阵行免疫细胞为标准进行scale正态分布,scale后去除na空白值
exprSet.map1=t(scale(t(exp_result1)))%>%
  na.omit()

#进行scale后热图发现还是有极值出现，遂设定最大最小值情况
table(abs(exprSet.map1)>2)
exprSet.map1[exprSet.map1>=2]=2
exprSet.map1[exprSet.map1<=-2]=-2


#设定分组标注，自定义分组标注颜色
annotation_col1=data.frame(CellType=c(rep("normal",10),rep("RA",13)))
row.names(annotation_col1)=colnames(exp_result1)
col=list(CellType=c(normal = "#87CEFA", RA = "#FFB6C1"))

text_rows=sample(rownames(exp_result1))

library(ggheatmap)
setwd("D:\\bioarticle/6.22/figure/CIBERSORT immun infiltration/")
pdf("GSE55457_ggheatmap.pdf",width = 15,height = 10)
ggheatmap=ggheatmap(exprSet.map1,cluster_rows = T,cluster_cols =F,##不设置聚类树，分组就不会被隔开
                    color = colorRampPalette(c("#00CED1",
                                               "#AFEEEE",
                                               "#FFFFF0",
                                               "#FFA07A",
                                               "#FF4500"))(2000),
                    annotation_cols = annotation_col1,
                    annotation_color = col,
                    text_show_rows = text_rows,
                    text_show_cols = F,
                    cluster_num = c(4),
                    tree_color_rows = c("#FF69B4",
                                        "#ADFF2F","#87CEEB",
                                        "#9370DB"))

ggheatmap%>%
  ggheatmap_theme(1:2,
                  theme = list(
                    theme(axis.text.y = element_text(color = "black",face = "bold",size = 15)),
                    theme(legend.title = element_text(face = "bold",size = 20))
                  ))
dev.off()



####22种免疫细胞组间比较，进行差异分析并箱线图可视化
##1.构建差异分析矩阵
group_list=c(rep("normal",10),rep("RA",13))
group_list=factor(group_list)%>%
  relevel(ref = "normal")

design=model.matrix(~group_list)
colnames(design)=levels(group_list)
row.names(design)=colnames(exp_result1)

##2.线性拟合模型构建，以及贝叶斯平滑处理
fit=lmFit(exp_result1,design)%>%
  eBayes()

##3.差异结果输出
allDiff_RME=topTable(fit,
                     coef = 2,
                     adjust.method = "fdr",
                     number = Inf)

##4.筛选具有显著性差异的细胞，p<0.05。
cell_p=subset(allDiff_RME,allDiff_RME$P.Val<0.05)


##5.提取差异细胞矩阵信息，并粗处理
RME_data=exp_result1[row.names(cell_p),]%>%
  t()%>%
  as.data.frame()

RME_data$group <- group_list
RME_data$sample <- row.names(RME_data)


##6.融合数据
library(reshape2)
RME_New=melt(RME_data)

colnames(RME_New)=c("Group","Sample","Celltype","Composition")  #设置行名


##7.画图
??stat_compare_means
library(ggpubr)
setwd("D:\\bioarticle/6.22/figure/CIBERSORT immun infiltration/")
pdf("GSE55457_boxplot.pdf",width = 10,height = 5)
ggplot(RME_New, aes(x = Celltype, y = Composition))+ 
  labs(y="Cell composition",x= NULL)+  
  geom_boxplot(aes(color = Group),##箱图线条色与填充色以分组为标准
               position=position_dodge(0.5),
               width=0.5,##误差线（上下横线长度）
               alpha = 0)+  ####设置箱子填充色透明度，0为完全透明，1为完全不透明(因为没有设置填充色标准，此处也可以不设置)
  scale_color_manual(values = c("#87CEFA", "#FFB6C1"))+##手动设置线条颜色（包括离散点）
  theme_classic() + 
  stat_compare_means(aes(group =  Group),##设置P值标志
                     label = "p.signif",
                     method = "wilcox.test",
                     hide.ns = T)+
  theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
        axis.title = element_text(size = 12,color ="black"), ##设置坐标轴标题字体及尺寸
        axis.text = element_text(size= 12,color = "black"),##设置坐标轴字体及尺寸
        panel.grid.minor.y = element_blank(),##用于设置网格线
        panel.grid.minor.x = element_blank(),##用于设置网格线
        axis.text.x = element_text(angle = 45, hjust = 1 ),##x轴文字倾斜角度
        panel.grid=element_blank(),
        legend.position = "top",##图例所在位置
        legend.text = element_text(size= 12),##修改图例文字大小
        legend.title= element_text(size= 12))+##修改图列标题大小
  geom_jitter(aes(color=Group),shape=16,
              stat = "identity",
              position = position_jitterdodge(0.2),
              size=2)##加离散点，以分组标准设置颜色
dev.off()

