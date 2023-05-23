Sys.setenv(LANGUAGE="en")
options(stringsAsFactors = F)
rm(list = ls())

###读取预先处理过的GSE55235数据集（raw原始数据）的
###所需矩阵数据
setwd("D:\\bioarticle/6.22/Analysis of immune infiltration in rheumatoid arthritis/CIBERSORT analyse/GSE55235/")
load("CIBERSORT_need.Rdata")
exp1=exp.pl1.6
colnames(exp1)=colnames(exp1.6)

###矩阵处理
###处理要求：表达矩阵格式第一列是基因名，第一行是样品名，
###不能有重复基因名，第一列列名不能有空白。
###矩阵中不能存在空白或NA值，不要对表达量取log2
##之前处理的是原始数据raw，并且已经用rma（）标准化过了
##注：rma使用pm信号，会对exp数据进行log2处理，
##所以此处需要对log2处理过的矩阵进行2^返回处理
setwd("D:\\bioarticle/6.22/Analysis of immune infiltration in rheumatoid arthritis/CIBERSORT analyse/GSE55235/three file/")
exp2=2^exp1
fwrite(exp2,file = "DATA.txt",sep = "\t",
       row.names = T,col.names = T)
DATA=fread("DATA.txt",header = T,data.table = F,
           sep = "\t")###打开保存的exp2矩阵看是否符合矩阵处理要求，发现不符合则重新处理
library(tidyverse)
DATA=rename(DATA,"Gene Symbol"="V1")
fwrite(DATA,file = "DATA.txt",sep = "\t",
       row.names = F,col.names = T)###处理完原exp2数据再次保存为原文件名文件

lm22=fread("LM22.txt",header = T,data.table = F,
           sep = "\t")###打开从nature文章下载的参考数据集（该参考数据集，提供22种免疫细胞亚型的基因表达特征集：LM22。
                      ###LM22.txt获取方法：
                      ###在Cibersort论文中下载Supplementry table 1
                      ###https://www.nature.com/articles/nmeth.3337#MOESM207)

##将LM22处理成满足上述矩阵处理要求的矩阵，并保存为原文件名文件
lm22.1=lm22[,-c(1:2)]
lm22.2=lm22.1%>%
  distinct(`Gene symbol`,.keep_all = T)%>%
  na.omit()
fwrite(lm22.2,file = "LM22.txt",sep = "\t",
       row.names = F,col.names = T)


###在R中新建R Script，复制以下网址中代码，保存为“Cibersort.R”
###https://rdrr.io/github/singha53/amritr/src/R/supportFunc_cibersort.R

###将”DATA.txt",“LM22.txt","Cibersort.R"保存在同一文件夹，运行Cibersort的代码如下：
setwd("D:\\bioarticle/6.22/Analysis of immune infiltration in rheumatoid arthritis/CIBERSORT analyse/GSE55235/three file/")
source("Cibersort.R")
result=CIBERSORT("LM22.txt","DATA.txt",perm = 1000,QN=T)
#perm置换次数=1000，QN分位数归一化=TRUE                
#QN如果是芯片设置为T，如果是测序就设置为F  

##运行时可能会跳出：Error in (function (cond) : 
## error in evaluating the argument 'x' in selecting a method for function 'sort': there is no package called ‘e1071’
##安装“e1071"R包即可，安装后再次返回此处运行result
BiocManager::install("e1071")
a

###注意：Cibersort结果的默认文件名为CIBERSORT-Results.txt，
###在同一文件夹下进行第二次运算会覆盖第一次得到的文件，
###建议在每一次运算之后对文件重命名。


####此步骤为处理ggplot热图所需数据格式（宽数据变长数据），
####但后来发现该方法不适合数据内容和图片格式要求，
####但所处理数据和步骤仍给予保留，以供将来参考
###处理result矩阵信息，把宽数据变为长数据，并进行热图展示
result=as.data.frame(result)
result1=rownames_to_column(result,"sample")
result2=result1[,1:23]
result3=pivot_longer(result2,-c("sample"),
                     names_to = "Type",
                     values_to = "Composition")

install.packages("scico")
library(scico)
library(ggplot2)
p1=ggplot(result3,aes(sample,Type))+
  geom_tile(aes(fill=Composition),colour="white",size=0.5)+
  theme_minimal()+
  scico::scale_fill_scico(palette = "bilbao")
####============分割


####ggheatmap免疫浸润热图展示
###发现ggplot画热图不行，遂用ggheatmap画热图，
###并预先处理得到所需数据格式（矩阵）
exp_result1=as.matrix(t(result[,1:22]))

##以矩阵行免疫细胞为标准进行scale正态分布,scale后去除na空白值
exprSet.map1=t(scale(t(exp_result1)))%>%
  na.omit()

#进行scale后热图发现还是有极值出现，遂设定最大最小值情况
table(abs(exprSet.map1)>1)
exprSet.map1[exprSet.map1>=1]=1
exprSet.map1[exprSet.map1<=-1]=-1

#设定分组标注，自定义分组标注颜色
annotation_col1=data.frame(CellType=c(rep("normal",10),rep("RA",10)))##忘记读取临床信息，就回看DEGs analyse的脚本复制过来这行代码，复制后记得修改
row.names(annotation_col1)=colnames(exp_result1)
col=list(CellType=c(normal = "#87CEFA", RA = "#FFB6C1"))

text_rows=sample(rownames(exp_result1))

BiocManager::install("ggheatmap")
a
library(ggheatmap)
setwd("D:\\bioarticle/6.22/figure/CIBERSORT immun infiltration/")
pdf("GSE55235_ggheatmap.pdf",width = 15,height = 10)
ggheatmap=ggheatmap(exprSet.map1,cluster_rows = T,cluster_cols = F,
          color = colorRampPalette(c("#00CED1",
                                     "#AFEEEE",
                                     "#FFFFF0",
                                     "#FFA07A",
                                     "#FF4500"))(2000),
          annotation_cols = annotation_col1,
          annotation_color = col,
          text_show_rows = text_rows,
          text_show_cols = F,
          cluster_num = c(4,5),
          tree_color_cols = c("#FF0000","#32CD32",
                              "#00BFFF","#8A2BE2"),
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
?pdf
?ggheatmap_theme


####22种免疫细胞组间比较，进行差异分析并箱线图可视化
##1.构建差异分析矩阵
group_list=c(rep("normal",10),rep("RA",10))
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
?geom_boxplot

?stat_compare_means
library(ggpubr)
setwd("D:\\bioarticle/6.22/figure/CIBERSORT immun infiltration/")
pdf("GSE55235_boxplot.pdf",width = 10,height = 5)
ggplot(RME_New, aes(x = Celltype, y = Composition))+ 
  labs(y="Cell composition",x= NULL)+  
  geom_boxplot(aes(color=Group),##箱图线条色与填充色以分组为标准
               position=position_dodge(0.5),
               width=0.5,##误差线（上下横线长度）
               alpha=0)+ ####设置箱子填充色透明度，0为完全透明，1为完全不透明(因为没有设置填充色标准，此处也可以不设置)
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


#引用R包参考文献
install.packages("pacman")
library(pacman)
p_cite("ggheatmap")
p_cite("ggpubr")
