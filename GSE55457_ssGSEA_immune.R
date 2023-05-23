Sys.setenv(LANGUAGE="en")
options(stringsAsFactors = F)
rm(list = ls())

###ssGSEA是通过R包GSVA去实现的，在分析之前，我们需要安装或加载GSVA包，
library(GSVA)

###1.读入数据并粗处理：
setwd("D:\\bioarticle/6.22/Analysis of immune infiltration/ssGSEA/GSE55457/")
load("ssGSEA_need.rdata")
exp=exp.pl.6

#读入文献当中定义的28种免疫细胞的基因集合
install.packages("read.xlsx")
library(readxl)
marker=read_xlsx("metagene list.xlsx",col_names = T,sheet = 1)[,1:2]
?read_xlsx
colnames(marker)=c("Metagene","Cell type")
marker=marker[-c(1:2),]

###2.构建免疫细胞Marker基因集
geneset=split(marker$Metagene,marker$`Cell type`)

###3.ssGSEA:
res=gsva(as.matrix(exp),
         geneset,
         method="ssgsea",
         kcdf="Gaussian",
         mx.diff=F,
         verbose=F)
##4.Min-Max标准化：此步结果resm暂不用,仅为将来技术参考
resm=res
for (i in colnames(res)) {
  resm[,i]=(res[,i]-min(res[,i]))/(max(res[,i])-min(res[,i]))
}

###5.提取类风湿关节炎GSE55457关节炎组样本ssGSEA结果矩阵信息
ra_res=res[,11:23]

###6.进行热图并初步层次聚类，然后简单暴力分组
library(pheatmap)
p=pheatmap(ra_res,scale = "row")
hc=cutree(p$tree_col,7)
ac=as.data.frame(as.character(hc))
rownames(ac)=colnames(ra_res)
pheatmap(ra_res,annotation_col = ac,
         scale = "row")

###7.把暴力分组调整为合理的免疫细胞高低分组,高表达为Cluster1，低表达为Cluster2
group_list=ifelse(hc<3,"Cluster2","Cluster1")
table(group_list)
ac=as.data.frame(group_list)
row.names(ac)=colnames(ra_res)
colnames(ac)="Cluster"

###8.热图最后绘图及参数调整
##用之前的数据尝试下画图，发现分组被割开
##为防止分组被聚类隔开，遂用自定义聚类分支顺序，按照分组人为指定热图聚类顺序排序
Cluster2_num=c(1,2,3,4,5,11)
Cluster1_num=c(6,7,8,9,10,12,13)
ra_res=ra_res[,c(Cluster1_num,Cluster2_num)]


?pheatmap
cluster_hot=pheatmap(ra_res,annotation_col = ac,
                     scale = "row",
                     show_colnames = T,
                     annotation_colors = list(
                       Cluster=c(Cluster1="#FF6347",Cluster2="#00CED1")),
                     color = colorRampPalette(c("#1E90FF","#87CEFA","#E1FFFF",
                                                "#FFFFF0","#FFDEAD",
                                                "#FF4500", "#FF0000"))(1000),
                     border_color = c("#DCDCDC"),
                     cellwidth = 40,cellheight = 20,
                     main = "GSE55457",
                     fontsize = 20,
                     fontsize_row = 15,
                     fontsize_col = 10,
                     cluster_cols = F)

#pheatmap热图需要自定义函数保存pdf
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
setwd("D:\\bioarticle/6.22/figure/ssGSEA immun infiltration/")

save_pheatmap_pdf(cluster_hot,"GSE55457_cluster_pheatmap.pdf",width = 15,height = 11)



#####DERATGs在亚型间的差异分析
###1.提取所含DERATG目标基因的矩阵信息
setwd("D:\\bioarticle/6.22/")
library(data.table)
DERATGs=fread("DERATGs1.csv",header = F,data.table = F)
exp_deratgs=exp[DERATGs$V1,]

###2.对DERATGs目标基因表达矩阵进行一次提取，按照cluster分组顺序进行样本排列
group_clust1=subset(ac,ac$Cluster=="Cluster1")
group_clust2=subset(ac,ac$Cluster=="Cluster2")

exp_deratgs=exp_deratgs[,c(rownames(group_clust1),
                           rownames(group_clust2))]


###3.构建差异分析矩阵，进行差异分析
library(tidyverse)
group_clust_list=c(rep("Cluster1",7),rep("Cluster2",6))%>%
  factor()%>%
  relevel(ref = "Cluster1")

design_clust=model.matrix(~group_clust_list)
colnames(design_clust)=levels(group_clust_list)
row.names(design_clust)=colnames(exp_deratgs)

library(limma)
fit=lmFit(exp_deratgs,design_clust)%>%
  eBayes()
allDiff=topTable(fit,coef = 2,adjust.method = "fdr",number = Inf)

###4.选取具有显著性差异的DERATGs，p<0.05
p_allDiff=subset(allDiff,allDiff$P.Value<0.05)

###5.提取显著性差异目标基因的矩阵信息，并粗处理
RME_data=exp_deratgs[row.names(p_allDiff),]%>%
  t()%>%
  as.data.frame()

RME_data$group <- group_clust_list
RME_data$sample <- row.names(RME_data)

##6.融合数据
library(reshape2)
RME_New=melt(RME_data)

colnames(RME_New)=c("Subtype","Sample","Genes","Expression")  #设置行名


##7.画图
library(ggpubr)
library(ggplot2)
setwd("D:\\bioarticle/6.22/figure/ssGSEA immun infiltration/")
pdf("GSE55457_ssGSEA_boxplot.pdf",width = 20,height = 10)
ggplot(RME_New, aes(x=Genes,y=Expression))+ 
  labs(y="Expression",x= NULL)+  
  geom_boxplot(aes(color=Subtype,fill=Subtype),##箱图线条色与填充色以分组为标准
               size=1,##误差线粗细
               width=1,##误差线（上下横线长度）
               alpha=0,##设置箱子填充色透明度，0为完全透明，1为完全不透明
               notch=F)+ ##是否为箱子添加槽口
  scale_color_manual(values = c("#FF6347","#00CED1"))+##手动调色
  theme_classic() + 
  theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
        axis.title = element_text(size = 30,color ="black"), ##设置坐标轴标题字体及尺寸
        axis.text = element_text(size= 30,color = "black"),##设置坐标轴字体及尺寸
        axis.ticks.length = unit(10.5,"pt"),##坐标轴刻度线朝向及长度，pt是单位
        axis.ticks = element_line(size = 1),##坐标轴刻度宽度
        axis.line = element_line(size = 1),##坐标轴线宽度
        panel.grid.minor.y = element_blank(),##用于设置网格线
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1 ),##x轴文字倾斜角度
        panel.grid=element_blank(),
        legend.position = "top",##图例所在位置
        legend.text = element_text(size= 20),##修改图例文字大小
        legend.title= element_text(size= 30),##修改图列标题大小
        legend.key.size = unit(30,"pt"))+##修改图例尺寸
  geom_jitter(aes(color=Subtype),shape=16,
              stat = "identity",
              position = position_jitterdodge(0.2),
              size=5)##加离散点，以分组标准设置颜色

dev.off()

setwd("D:\\bioarticle/6.22/Analysis of immune infiltration/ssGSEA/GSE55457/")
