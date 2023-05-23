Sys.setenv(LANGUAGE="en")
options(stringsAsFactors = F)
rm(list=ls())

####GSE55457 GSEA及GSVA分析
##加载基因集矩阵、差异分析及临床表型数据
setwd("D:\\bioarticle/6.22/analyse data and script/GSEA and GSVA analyse data/GSE55457/")
load("GSE55457_GSEA_GSVA_need.rdata")


####GSEA
##加载GSEA分析所需R包
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(GSEABase)
keytypes(org.Hs.eg.db)

##添加entrezid列
#symbol转entrez ID：
symbol=row.names(allDiff2)
entrez=bitr(symbol,
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = org.Hs.eg.db)

##基因排序
DEG=merge(entrez,allDiff2,by.x="SYMBOL",by.y=0)
colnames(DEG)
library(tidyverse)
data_all_sort=DEG%>%
  arrange(desc(logFC))
gene_list=data_all_sort$logFC
names(gene_list)=as.character(data_all_sort$ENTREZID)
head(gene_list)

##得出GSEA富集结果
egmt=gseKEGG(geneList = gene_list,
             organism = "hsa",
             minGSSize = 10,
             maxGSSize = 500,
             pvalueCutoff = 0.05,
             pAdjustMethod = "BH",
             verbose = FALSE,
             eps = 0)

#提取富集结果
egmt.result=egmt@result

#以FDR<0.25且p<0.05作为标准，筛选显著富集的通路
egmt.result2=subset(egmt.result,
                    egmt.result$qvalues<0.25&egmt.result$pvalue<0.05)


##按照enrichment score从高到低排序，便于查看富集通路
egmt.sort=egmt.result2[order(egmt.result2$enrichmentScore,decreasing = T),]


##画ES排前10的高表达通路
setwd("D:\\bioarticle/6.22/figure/GSEA GSVA_enrich/GSE55457/")
pdf("GSE55457_up_kegg.pdf",width = 13,height = 10)
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
pdf("gse55457_down_kegg.pdf",width = 13,height = 10)
gseaplot2(egmt,
          geneSetID = row.names(egmt.sort)[46:55],
          subplots = 1:2,
          color =c("#6495ED","#00BFFF","#2F4F4F"
                   ,"#40E0D0","#006400","#191970",
                   "#4169E1","#00CED1","#3CB371",
                   "#7CFC00"),
          base_size = 23)
dev.off()

###制作GSE55457-GSEA分析结果（可视化部分）表格基本信息table
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
gsea_result2=gsea_result1[c(1:10,46:55),c("ID","Description",
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
          caption="Table 6.GSEA-KEGG analysis results of GSE55457 (top 10 terms of up-regulated and down-regulated pathways were listed)",
          align="left",
          tfoot="GSEA gene set enrichment analysis
                 KEGG kyoto encyclopedia of genes and genomes",
          rgroup = c("Up-regulated","Down-regulated"),
          n.rgroup = c(10,10))
#========================分割线========================
##决定还是把GSE55235和GSE55457的gsea、gsva数据合并再画table
library(tidyverse)
gsea_result3=mutate(gsea_result2,"GEO accession"="GSE55457",.before =1 )

#保存数据
setwd("D:\\bioarticle/6.22/analyse data and script/GSEA and GSVA analyse data/")
write.csv(gsea_result3,"gse55457_GSEA_table.csv")




###GSVA
##导入分析参考基因集
getwd()
setwd("D:\\bioarticle/6.22/")
kegg.symbol.gmt=read.gmt("c2.cp.kegg.v7.5.1.symbols (1).gmt")
##切割参考基因集对象，为后续GSVA分析做准备
symbol.gmt=split(kegg.symbol.gmt$gene,kegg.symbol.gmt$term)


library(GSVA)
##将基因样本表达矩阵从数据框格式转化为矩阵格式，为gsva分析做准备
expr=as.matrix(exp.pl.6)
##gsva分析
gsva=gsva(expr,symbol.gmt,kcdf="Gaussian",method="gsva")



##根据GEO样本名分组建立差异分析矩阵，进行差异分析
group_list=c(rep("normal",10),rep("RA",13))
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
setwd("D:\\bioarticle/6.22/analyse data and script/GSEA and GSVA analyse data/GSE55457/")
write.csv(gsva.allDiff,file = "GSE55457_gsva_allDiff.csv")

#筛选出高低表达的差异通路
dif=gsva.allDiff
dif$significance="stable"
log2(1.2)
dif$significance[dif$logFC>=0.263&dif$adj.P.Val<0.05]="up"
dif$significance[dif$logFC<=-0.263&dif$adj.P.Val<0.05]="down"
table(dif$significance)     
#运行后结果 
down stable     up 
1    171     14 
#筛取down和up的通路
dif.up=subset(dif,dif$significance=="up")
dif.down=subset(dif,dif$significance=="down")
dif.subset=rbind(dif.up,dif.down)

pheatmap.dif=gsva[row.names(dif.subset),]


###画热图
annotation_col=data.frame(Tissuetype=c(rep("normal",10),
                                       rep("RA",13)))
row.names(annotation_col)=colnames(exp.pl.6)

##以矩阵行通路为标准进行scale正态分布
exprSet.map1=t(scale(t(pheatmap.dif)))


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
getwd()
setwd("D:\\bioarticle/6.22/figure/GSEA GSVA_enrich/GSE55457/")


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

save_pheatmap_pdf(hot,"GSE55457_pheatmap.pdf",width = 20,height = 30)


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
          caption="Table 8.GSVA difference analysis results of GSE55457 (all terms of up-regulated and down-regulated pathways were listed)",
          align="left",
          tfoot="GSVA gene set variation analysis ",
          rgroup = c("Up-regulated","Down-regulated"),
          n.rgroup = c(14,1))
#=====================分割线===================================
library(tidyverse)
gsva_diff_table1=mutate(gsva_diff_table,
                        "GEO accession"="GSE55457",
                        "Description"=rownames(gsva_diff_table),
                        .before = 1)
#保存数据
setwd("D:\\bioarticle/6.22/analyse data and script/GSEA and GSVA analyse data/")
write.csv(gsva_diff_table1,"gse55457_GSVA_table.csv")
