Sys.setenv("LANGUAGE"="en")
options(stringsAsFactors = F)
rm(list = ls())

setwd("D:\\bioarticle/6.22/analyse data and script/GSEA and GSVA analyse data/")
##画GSEA的table
library(data.table)
gsea_55235=fread("gse55235_GSEA_table.csv",header = T,data.table = F)
gsea_55235=gsea_55235[,-1]

gsea_55457=fread("gse55457_GSEA_table.csv",header = T,data.table = F)
gsea_55457=gsea_55457[,-1]

gsea.table=rbind(gsea_55235,gsea_55457)
gsea.table=gsea.table[,-7]

library(pixiedust)
library(tidyverse)
library(kableExtra)

dust(gsea.table)%>%
  sprinkle_caption("Table 5.GSEA enrichment summary")%>%
  sprinkle_caption_number(caption_number = getOption("Table 5.GSEA enrichment summary",FALSE))%>%
  sprinkle_fn(cols = 5:6,fn = quote(pvalString(value)))%>%
  sprinkle(cols=4,round=3)%>%
  sprinkle_merge(cols = 1,merge = TRUE,rows = 1:20)%>%
  sprinkle_merge(cols = 1,merge = TRUE,rows = 21:40)%>%
  sprinkle_print_method(print_method = c("html"))

?sprinkle_border

##画GSVA的table
gsva.55235=fread("gse55235_GSVA_table.csv",header = T,data.table = F)
gsva.55235=gsva.55235[,-1]

gsva.55457=fread("gse55457_GSVA_table.csv",header = T,data.table = F)
gsva.55457=gsva.55457[,-1]

gsva.table=rbind(gsva.55235,gsva.55457)

library(pixiedust)
library(tidyverse)
library(kableExtra)

dust(gsva.table)%>%
  sprinkle_caption("Table 6.GSVA enrichment summary")%>%
  sprinkle_caption_number(caption_number = getOption("Table 6.GSVA enrichment summary",FALSE))%>%
  sprinkle_fn(cols = 5:6,fn = quote(pvalString(value)))%>%
  sprinkle(cols=3:4,round=3)%>%
  sprinkle_merge(cols = 1,merge = TRUE,rows = 1:39)%>%
  sprinkle_merge(cols = 1,merge = TRUE,rows = 40:54)%>%
  sprinkle_print_method(print_method = c("html"))
