Sys.setenv(LANGUAGE="en")
options(stringsAsFactors = F)
rm(list = ls())
setwd("D:\\bioarticle/6.22/")

###韦恩图.筛查靶基因
##1.数据准备
library(data.table)
DEGs1=fread("DEGs1.csv",header = T,data.table = F)
DEGs2=fread("DEGs2.csv",header = T,data.table = F)
RATGs=fread("GeneCards-SearchResults.csv",header = T,data.table = F)

x=list(GSE55235=DEGs1$V1,
       GSE55457=DEGs2$V1,
       Target=RATGs$`Gene Symbol`)
##2.画图并保存
install.packages("VennDiagram")
library(VennDiagram)
##2022-07-26修图运行时提示：载入需要的程辑包：grid
##载入需要的程辑包：futile.logger
BiocManager::install("grid")
a
library(grid)
BiocManager::install("futile.logger")
a
library(futile.logger)
?venn.diagram
setwd("D:\\bioarticle/6.22/figure/venn/")
venn.diagram(x,fill=c("violet","red","green"),
             alpha=0.5,
             main = "DEG_Target",
             main.fontface = "plain",
             main.fontfamily = "serif",
             sub.fontface = "plain",
             sub.fontfamily = "serif",
             main.cex = 2,
             filename = "Venn.tiff")
##venn.diagram函数只能输出png和tiff格式



##尝试UpSetR函数画韦恩图
install.packages("UpSetR") 
library(UpSetR)
UpSetR::upset(fromList(x),
              nsets = 3,order.by = "freq",
              point.size = 5,
              line.size = 1.3,
              mainbar.y.label = "DEG_Target",
              sets.x.label = "",
              mb.ratio = c(0.60,0.40),
              text.scale = c(2,2,2,2,2,2,2))
##呵呵，没画出来，等有时间再研究这个怎么画韦恩图



##3.取出重叠区靶基因
area=calculate.overlap(x)
typeof(area)
summary(area)
DERATGs=area$a5
setwd("D:\\bioarticle/6.22/")
write.csv(DERATGs,file = "DERATGs.csv")
write.table(DERATGs,file = "DGRATGs.txt")
load("deratgs.rdata.RData")
?summary
class(area)
library(data.table)
DERATGs1=fread("DERATGs1.csv",header = F,data.table = F)

###GO与KEGG富集分析clusterProfiler
BiocManager::install("clusterProfiler")
a
BiocManager::install("enrichplot")
a
BiocManager::install("GOplot")
a
library(clusterProfiler)
install.packages("viridis")
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

library(GOplot)
library(tidyverse)

Sys.setenv(LANGUAGE="en")
options(stringsAsFactors = F)
rm(list = ls())
setwd("D:\\bioarticle/6.22/")
load("enrich analyse.RData")


###数据准备，将symbol转换为entrezid

degs1=DEGs1
row.names(degs1)=degs1$V1
diff1=degs1[DERATGs,]

diff1$SYMBOL=diff1$V1
diff2=diff1[,-1]
keytypes(org.Hs.eg.db)
?bitr
gene.df=bitr(diff2$SYMBOL,
             fromType = "SYMBOL",
             toType = c("ENTREZID"),
             OrgDb = org.Hs.eg.db)

#太难了，不知道这是啥Enrich=gene.df%>%inner_join(diff2,by="SYMBOL")
class(Enrich)
rm(Enrich)
##GO富集，得出GO富集结果ego
ego=enrichGO(gene = gene.df$ENTREZID,##ENTREZID
             keyType = "ENTREZID",##输入基因类型
             OrgDb = org.Hs.eg.db,##导入背景基因
             ont = "all",##GO的种类，BP、CC、MF
             pAdjustMethod = "BH",##矫正的类型
             pvalueCutoff = 0.05,##p值过滤值
             readable = TRUE)##readable=TRUE或setReadable函数可以将ENTREZID转换回symbol
?enrichGO
head(ego)


###GO注释柱状图
setwd("D:\\bioarticle/6.22/figure/GO KEGG DO_enrich/")
pdf("key_GO_bar.pdf",width = 15,height=10)
barplot(ego,
        drop=TRUE,
        showCategory = 10,
        split="ONTOLOGY")+
  scale_y_discrete(labels=function(x)str_wrap(x,width = 80))+
  facet_grid(ONTOLOGY~.,scale ="free" )
?barplot

dev.off()

###GO气泡图
?dotplot
pdf("key_GO_dot.PDF",width = 15,height = 10)
dotplot(ego,
        showCategory=10,
        split="ONTOLOGY",
        label_format=100)+
  facet_grid(ONTOLOGY~.,scale="free")

dev.off()
setwd("D:\\bioarticle/6.22/figure/enrich/")
setwd("D:\\bioarticle/6.22/")

###GO网络图gene-concept network
?cnetplot
##1.获取gene_list
gene_list=diff2$logFC
##2.命名
names(gene_list)=gene.df$ENTREZID
##3.排序很重要
gene_list=sort(gene_list,decreasing = TRUE)
View(gene_list)
##4.画图
setwd("D:\\bioarticle/6.22/figure/GO KEGG DO_enrich/")
pdf("key_GO_cnetplot1.pdf",width = 15,height = 10)
cnetplot(ego,foldChange = gene_list,circular=TRUE,
         colorEdge=TRUE,showCategory = 10)
BiocManager::install("ggnewscale")
a
library(ggnewscale)
dev.off()


###制作GO分析结果（可视化部分）表格基本信息table
#数据准备
ego_result=ego@result
colnames(ego_result)
ego_result1=ego_result[,c("ONTOLOGY","ID","Description",
                          "p.adjust","geneID","Count")]
ego_result2=rbind(filter(ego_result1,ONTOLOGY=="BP")[1:10,],
                  filter(ego_result1,ONTOLOGY=="CC")[1:10,],
                  filter(ego_result1,ONTOLOGY=="MF")[1:10,])%>%
  na.omit()
setwd("D:\\bioarticle/6.22/analyse data and script/Ven and GO KEGG DO enrich analyse data/")
write.table(ego_result2,file = "GO analysis results of DERATGs",
            row.names = FALSE,sep = ",")


str(ego_result2)#查看分类变量数据类型
#分类变量是字符型，接下来进行变量因子转换
ego_result2$ONTOLOGY=factor(ego_result2$ONTOLOGY,
                            ordered = T,
                            levels = c("BP","CC","MF"))

##==============================分割线===========================================
#安装加载R包
#制作表格
#方法一
library(Gmisc)
library(Rcpp)
library(htmlTable)
library(kableExtra)
install.packages("magick")
library(magick)
library(tidyverse)
install.packages("PhantomJS")
install_phantomjs(version = "2.1.1",
                  baseURL = "https://github.com/wch/webshot/releases/download/v0.3.1/")

?htmlTable

htmlTable(ego_result2,
          rname=F,#设置行名不显示
          caption="Table 2.GO analysis results of DERATGs(top 10 terms of BP，top 10 terms of MF and all terms of CC were listed)",
          align="left",
          tfoot="GO gene ontology, DERATGs differentially expressed rheumatoid arthritis target genes
          BP biological process,CC cellular component,MF molecular function")


###======================================================分割线===========================================
#方法二"kableExtra"
install.packages("kableExtra")
library(kableExtra)
install.packages("rmarkdown")
library(rmarkdown)
setwd("D:\\bioarticle/6.22/analyse data and script/Ven and GO KEGG DO enrich analyse data/")

kbl(ego_result2[,2:6], 
    row.names = F,
    caption = "Table 2.GO analysis results of DERATGs") %>%
  kable_paper("striped", full_width = F) %>%#设置页面样式及是否真实宽度
  pack_rows("BP", 1, 10) %>%
  pack_rows("CC", 11, 12)%>%
  pack_rows("CC", 13, 22)%>%
  footnote(general = "GO gene ontology, DERATGs differentially expressed rheumatoid arthritis target genes",
           title_format = "italic")%>%
  kable_styling(font_size = 21,#设置字体
                fixed_thead = F,#表格太长，滚动时将标题行固定在顶部
                position = "center")
###===============================分割线===============================================

#方法三
install.packages("pixiedust")
library(pixiedust)
library(tidyverse)
library(kableExtra)
#之前方法用的ego_result2中有因子变量，此处需要转化为原来的字符变量才行
str(ego_result2)
ego_result3=ego_result2
ego_result3$ONTOLOGY=as.character(ego_result3$ONTOLOGY)
setwd("D:\\bioarticle/6.22/analyse data and script/Ven and GO KEGG DO enrich analyse data/")
dust(ego_result3[,-5])%>%
  sprinkle_caption("Table 2.GO enrichment summary")%>%
  sprinkle_caption_number(caption_number = getOption("Table 2.GO enrichment summary",FALSE))%>%
  sprinkle_fn(cols = 4,fn = quote(pvalString(value)))%>%
  sprinkle_print_method(print_method = c("html"))


?sprinkle_print_method
?sprinkle_caption
?sprinkle_border_collapse.default
?sprinkle_caption_number
##=======================分割线====================================
##最后定了用方法三，数值看起来清爽一些
####输出到画板Viewer显示，然后通过网页打开，网页另存为html文件，再用word打开修改表格
####试了一天也没找到方法直接输出成好看的PDF格式文件，只能先word手动操作了，简直呵呵





###KEGG富集分析
?enrichKEGG
ekegg=enrichKEGG(gene = gene.df$ENTREZID,
                keyType = "kegg",
                organism = "hsa",
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH")


###kegg注释柱状图
setwd("D:\\bioarticle/6.22/figure/GO KEGG DO_enrich/")
pdf("key_kegg_bar.pdf",width = 15,height=10)
barplot(ekegg,
        drop=TRUE,
        showCategory = 7,
        label_format = 100)

dev.off()

###kegg气泡图
getwd()
pdf("key_kegg_dotplot.pdf",width = 15,height = 10)
dotplot(ekegg,
        showCategory=7,
        label_format=100)
dev.off()

###kegg网络图
##gene_list和之前一样
##kegg返回Gene Symbol名称
library(DOSE)
ekegg_rs=DOSE::setReadable(ekegg,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
pdf("key_kegg_cnetplot.pdf",width = 15,height = 10)
cnetplot(ekegg_rs,foldChange = gene_list,
         showCategory = 7,circular=TRUE,
         colorEdge=TRUE)
dev.off()

##kegg通路图，得到pi3k,需要取浏览器输入网址可视化
kegg1=ekegg@result
pi3k=browseKEGG(kegg1,"hsa04151")


###制作KEGG分析结果（可视化部分）表格基本信息table
#数据准备
colnames(kegg1)#kegg1用的是entrezid，我们需要gene symbol的内容
kegg_result=ekegg_rs@result
kegg_result1=kegg_result[1:7,c("ID","Description","p.adjust","qvalue","geneID","Count")]
##=============================分割线=================================
#此处备用，代码先留着
#加载R包，画table
library(htmlTable)
htmlTable(kegg_result1,
          rname=F,#设置行名不显示
          caption="Table 3.KEGG analysis results of DERATGs (top 7 terms were listed)",
          align="left",
          tfoot="KEGG kyoto encyclopedia of genes and genomes
          DERATGs differentially expressed rheumatoid arthritis target genes")
##=======================分割线================================
install.packages("pixiedust")
library(pixiedust)
library(tidyverse)
library(kableExtra)

dust(kegg_result1[,-5])%>%
  sprinkle_caption("Table 3. KEGG enrichment summary")%>%
  sprinkle_caption_number(caption_number = getOption("Table 3. KEGG enrichment summary",FALSE))%>%
  sprinkle_fn(cols = 3:4,fn = quote(pvalString(value)))%>%
  sprinkle_print_method(print_method = c("html"))



###Disease Ontology富集分析及可视化
BiocManager::install("DOSE")
a
library(DOSE)
edo=enrichDO(gene = gene.df$ENTREZID,
             ont = "DO",
             pvalueCutoff = 0.05,
             pAdjustMethod = "BH",
             qvalueCutoff = 0.05,
             readable = TRUE)
###DO注释柱状图
getwd()
?barplot
pdf("key_DO_bar.pdf",width = 15,height = 10)
edo1=edo@result
barplot(edo,showCategory = 10,
        label_format = 100)
dev.off()
setwd("D:\\bioarticle/6.22/figure/enrich/")
load("enrich analyse.RData")

###DO气泡图
pdf("key_DO_dot.pdf",width = 15,height = 10)
dotplot(edo,showCategory=10,
        label_format=100)
dev.off()

###DO网络图
##gene_list和之前一样
pdf("key_DO_cnetplot.pdf",width = 15,height = 10)
cnetplot(edo,foldChange = gene_list,
         showCategory = 10,circular=TRUE,
         colorEdge=TRUE)
dev.off()

setwd("D:\\bioarticle/6.22/")

###制作DO分析结果（可视化部分）表格基本信息table
#数据准备
colnames(edo1)
do_result=edo1[1:10,c("ID","Description","p.adjust","qvalue","geneID","Count")]

#============================分割线====================================================
#此处备用，代码先留着
#加载R包，画table
library(htmlTable)
htmlTable(do_result,
          rname=F,#设置行名不显示
          caption="Table 4.DO analysis results of DERATGs (top 10 terms were listed)",
          align="left",
          tfoot="DO disease ontology analysis
          DERATGs differentially expressed rheumatoid arthritis target genes")
#============================分割线===========================================

install.packages("pixiedust")
library(pixiedust)
library(tidyverse)
library(kableExtra)

dust(do_result[,-5])%>%
  sprinkle_caption("Table 4. DO enrichment summary")%>%
  sprinkle_caption_number(caption_number = getOption("Table 4. DO enrichment summary",FALSE))%>%
  sprinkle_fn(cols = 3:4,fn = quote(pvalString(value)))%>%
  sprinkle_print_method(print_method = c("html"))





















####============================分割线===================================
####下面的GSEA和GSVA只是尝试分析，具体见GSEA、GSVA的正式分析
####GSEA与GSVA分析及可视化
kegg.gmt=read.gmt("c2.cp.kegg.v7.5.1.symbols (1).gmt")
kegg_list=split(kegg.gmt$gene,kegg.gmt$term)
kegg.gmt.1=read.gmt("c2.cp.kegg.v7.5.1.entrez.gmt")

BiocManager::install("GSEABase")
a
library(GSEABase)
BiocManager::install("GSVA")
a
library(GSVA)
view(gene_list)

###GSEA分析
#先用DERATGs看一下
?GSEA
head(gene_list)
gsea.deratgs=GSEA(gene_list,
                  TERM2GENE = kegg.gmt.1,
                  pvalueCutoff = 1)
head(gsea.deratgs)
gsea.deratgs.1=gsea.deratgs@result

?gseaplot2
setwd("D:\\bioarticle/6.22/figure/GSEA GSVA_enrich/")
pdf("DERATGs.pdf",width = 15,height = 10)
gseaplot2(gsea.deratgs,
          geneSetID = gsea.deratgs.1$ID,
          subplots = 1:2,
          pvalue_table = T,
          base_size = 20)
dev.off()

###GSE55235的GSEA富集及可视化
####=================================分割线===================================




