Sys.setenv(LANGUAGE="en")
options(stringsAsFactors = F)
rm(list = ls())
setwd("D:\\bioarticle/6.22/Cytoscape data/")
library(data.table)

####string得出47个DERATGs所发生PPI互作关系的数目进行可视化

###cytoscape软件用的是string_interactions_short.tsv，
###此处PPI互作关系数目文件需要从string里重新下载
cyto_edge=fread("string_interactions.tsv",header = T,data.table = F)
?aggregate
cyto_edge.1=cyto_edge[,1:2]
colnames(cyto_edge.1)=c("node1","node2")
cyto_edge.1$edges=c(1)
cyto_edge.2=cyto_edge.1[,-2]

##处理后的data.frame文件进行aggregate去重复及edge数加和处理
edge1=aggregate(cyto_edge.2$edges,
                list(cyto_edge.2$node1),
                FUN = sum)

colnames(edge1)=c("node","edge")
#按照edge进行升序排序
edge2=edge1[order(edge1$edge,decreasing = F),]
row.names(edge2)=edge2$node

library(ggplot2)
#怎加渐变色对象，为后续柱状图赋色准备
colors=colorRampPalette(c("#4169E1","#87CEFA",#87CEFA",
                          "#7FFFAA","#FFA500",
                          "#FF0000","#C71585"))(47)

#画柱状图
setwd("D:\\bioarticle/6.22/figure/PPI/")
pdf("PPI_ggplot.pdf",width = 20,height = 15)
ggplot(edge2,aes(x=edge,y=node))+#确定坐标轴内容
  geom_bar(stat = "identity",color=colors,
           fill=colors)+##画初始图并上色
  labs(x="number of edges",y="name of nodes")+##设立坐标轴名称
  theme_classic()+##设定图的主题板式
  scale_y_discrete(limits=row.names(edge2))+#y轴排序
  geom_text(aes(label=edge),size=10)+#为柱状图添加计数标记
  theme(text = element_text(size = 30))#设置文本大小
dev.off()
?pdf
?ggplot

setwd("D:\\bioarticle/6.22/Cytoscape data/")
write.csv(edge2,file = "nodes.csv")
