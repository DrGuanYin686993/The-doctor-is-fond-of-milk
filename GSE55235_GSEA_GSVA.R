Sys.setenv(LANGUAGE="en")
options(stringsAsFactors = F)
rm(list = ls())

setwd("D:\\bioarticle/6.22/analyse data and script/GSEA and GSVA analyse data/GSE55235/")
getwd()
load("GSE55235_GSEA_GSVA_need.Rdata")
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(GSEABase)
keytypes(org.Hs.eg.db)

####GSE55235的GSEA与GSVA分析及可视化

###GSEA
##添加entrezid列
#symbol转entrez ID：
symbol=row.names(data1)
entrez=bitr(symbol,
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = org.Hs.eg.db)

##基因排序
DEG=merge(entrez,data1,by.x="SYMBOL",by.y=0)
colnames(DEG)
library(tidyverse)
data_all_sort=DEG%>%
  arrange(desc(logFC))
gene_list=data_all_sort$logFC
names(gene_list)=as.character(data_all_sort$ENTREZID)
head(gene_list)
?desc
?arrange

##得出富集结果
egmt=gseKEGG(geneList = gene_list,
             organism = "hsa",
             minGSSize = 10,
             maxGSSize = 500,
             pvalueCutoff = 0.05,
             pAdjustMethod = "BH",
             verbose = FALSE,
             eps = 0)
?gseKEGG

##提取富集结果
egmt.result=egmt@result

#以FDR<0.25且p<0.05作为标准，筛选显著富集的通路
egmt.result2=subset(egmt.result,
                    egmt.result$qvalues<0.25&egmt.result$pvalue<0.05)

##GSEA分析高表达通路
##按照enrichment score从高到低排序，便于查看富集通路
egmt.sort=egmt.result2[order(egmt.result2$enrichmentScore,decreasing = T),]


##画前10个上调通路
setwd("D:\\bioarticle/6.22/figure/GSEA GSVA_enrich/GSE55235/")
pdf("GSE55235_up_kegg.pdf",width = 13,height = 10)
gseaplot2(egmt,
          geneSetID = row.names(egmt.sort)[1:10],
          subplots = 1:2,
          color =c("#FF0000","#C71585","#FF4500"
                   ,"#FFD700","#D2691E","#800000",
                   "#D2691E","#CD853F","#FF6347",
                   "#CD5C5C"),
          base_size = 23)
dev.off()

##画后10个下调通路


pdf("gse55235_down_kegg.pdf",width = 13,height = 10)
gseaplot2(egmt,
          geneSetID = row.names(egmt.sort)[65:74],
          subplots = 1:2,
          color =c("#6495ED","#00BFFF","#2F4F4F"
                   ,"#40E0D0","#006400","#191970",
                   "#4169E1","#00CED1","#3CB371",
                   "#7CFC00"),
          base_size = 23)
dev.off()


###制作GSE55235-GSEA分析结果（可视化部分）表格基本信息table
##数据准备
#egmt返回Gene Symbol名称
library(DOSE)
library(org.Hs.eg.db)
egmt_symbol=DOSE::setReadable(egmt,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
gsea_result=egmt_symbol@result
#以FDR<0.25且p<0.05作为标准，筛选显著富集的通路
gsea_result1=filter(gsea_result,
                    gsea_result$qvalues<0.25&gsea_result$pvalue<0.05)
#按照enrichment score从高到低排序，便于查看富集通路
gsea_result1=gsea_result1[order(gsea_result1$enrichmentScore,decreasing = T),]
#提取top10上调通路和top10下调通路
colnames(gsea_result1)
gsea_result2=gsea_result1[c(1:10,65:74),c("ID","Description",
                                     "enrichmentScore",
                                     "pvalue","qvalues",
                                     "core_enrichment")]

#======================分割线===============================
#之后决定用pixiedust的函数画table，但此处代码仍保留
#加载R包，画table
library(htmlTable)
?htmlTable
htmlTable(gsea_result2,
          rname=F,#设置行名不显示
          caption="Table 5.GSEA-KEGG analysis results of GSE55235 (top 10 terms of up-regulated and down-regulated pathways were listed)",
          align="left",
          tfoot="GSEA gene set enrichment analysis
                 KEGG kyoto encyclopedia of genes and genomes",
          rgroup = c("Up-regulated","Down-regulated"),
          n.rgroup = c(10,10))
#========================分割线========================

install.packages("pixiedust")
library(pixiedust)
library(tidyverse)
library(kableExtra)

dust(gsea_result3)%>%
  sprinkle_caption("Table 5.GSEA enrichment summary")%>%
  sprinkle_caption_number(caption_number = getOption("Table 5.GSEA enrichment summary",FALSE))%>%
  sprinkle_fn(cols = 5:6,fn = quote(pvalString(value)))%>%
  sprinkle(cols=4,round=3)%>%
  sprinkle_merge(cols = 1,merge = TRUE,rows = 1:10)%>%
  sprinkle_merge(cols = 1,merge = TRUE,rows = 11:20)%>%
  sprinkle_print_method(print_method = c("html"))

??sprinkle_merge.default
?mutate
gsea_result3=mutate(gsea_result2,"GEO accession"="GSE55235",.before =1 )
#代码敲到这里，决定还是把GSE55235和GSE55457的gsea、gsva数据合并再画table
#保存数据
setwd("D:\\bioarticle/6.22/analyse data and script/GSEA and GSVA analyse data/")
write.csv(gsea_result3,"gse55235_GSEA_table.csv")





###GSVA
##加载预定基因集
setwd("D:\\bioarticle/6.22/")
kegg.symbol.gmt=read.gmt("c2.cp.kegg.v7.5.1.symbols (1).gmt")

#将kegg基因集切割，为后面GSVA分析做准备
symbol.gmt=split(kegg.symbol.gmt$gene,kegg.symbol.gmt$term)


colnames(exp.pl1.6)=colnames(exp1.6)
BiocManager::install("GSVA")
a
library(GSVA)
##要将表达矩阵转换为matrix格式
expr=as.matrix(exp.pl1.6)
##gsva分析
gsva=gsva(expr,symbol.gmt,kcdf="Gaussian",method="gsva")
##根据GEO样本名分组建立差异分析矩阵，进行差异分析
group_list=c(rep("normal",10),rep("RA",10))
group_list=factor(group_list)
group_list=relevel(group_list,ref = "normal")
design=model.matrix(~group_list)
colnames(design)=levels(group_list)
row.names(design)=colnames(gsva)

library(limma)
fit=lmFit(gsva,design)
fit=eBayes(fit)
gsva.allDiff=topTable(fit,coef = 2,adjust.method = "fdr",number = Inf)
#保存差异分析结果
setwd("D:\\bioarticle/6.22/analyse data and script/GSEA and GSVA analyse data/GSE55235/")
write.csv(gsva.allDiff,file = "GSE55235_gsva_allDiff.csv")
#筛选出高低表达的差异通路
dif=gsva.allDiff
dif$significance="stable"
log2(1.2)
dif$significance[dif$logFC>=0.263&dif$adj.P.Val<0.05]="up"
dif$significance[dif$logFC<=-0.263&dif$adj.P.Val<0.05]="down"
table(dif$significance)
dif.up=subset(dif,dif$significance=="up")
dif.down=subset(dif,dif$significance=="down")
dif.subset=rbind(dif.up,dif.down)

pheatmap.dif.subset=gsva[row.names(dif.subset),]


##画热图
annotation_col=data.frame(Tissuetype=c(rep("normal",10),
                                       rep("RA",10)))
row.names(annotation_col)=colnames(exp.pl1.6)

library(pheatmap)
getwd()
setwd("D:\\bioarticle/6.22/figure/GSEA GSVA_enrich/GSE55235/")

##以矩阵行基因为标准进行scale正态分布
exprSet.map1=t(scale(t(pheatmap.dif.subset)))
#进行scale后热图发现还是有极值出现，遂设定最大最小值情况
table(abs(exprSet.map1)>1.5)
exprSet.map1[exprSet.map1>=1.5]=1.5
exprSet.map1[exprSet.map1<=-1.5]=-1.5

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


?pheatmap
hot=pheatmap(exprSet.map1,
             annotation_col = annotation_col,
             annotation_colors = ann_colors,
             show_colnames = F,
             show_rownames = T,
             cluster_cols = col_cluster,
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

save_pheatmap_pdf(hot,"GSE55235_pheatmap.pdf",width = 20,height = 30)

getwd()

setwd("GSEA and GSVA analyse data/GSE55235/")



###制作GSVA富集数据差异分析结果（可视化部分）表格基本信息table
#数据准备
colnames(dif.subset)
gsva_diff_table=dif.subset[,c("logFC","AveExpr","P.Value","adj.P.Val")]

#======================分割线===============================
#之后决定用pixiedust的函数画table，但此处代码仍保留
#加载R包，画table
library(htmlTable)
?htmlTable
htmlTable(gsva_diff_table,
          caption="Table 7.GSVA difference analysis results of GSE55235 (all terms of up-regulated and down-regulated pathways were listed)",
          align="left",
          tfoot="GSVA gene set variation analysis ",
          rgroup = c("Up-regulated","Down-regulated"),
          n.rgroup = c(32,7))
#=====================分割线===================================

gsva_diff_table1=mutate(gsva_diff_table,
                        "GEO accession"="GSE55235",
                        "Description"=rownames(gsva_diff_table),
                        .before = 1)
#保存数据
setwd("D:\\bioarticle/6.22/analyse data and script/GSEA and GSVA analyse data/")
write.csv(gsva_diff_table1,"gse55235_GSVA_table.csv")
