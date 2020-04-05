load("Data/STEP2_de_result.RData")
load("Data/STEP1_input.RData")
source("utils.R")

#===模式生物的KEGG富集分析,enrichKEGG ===
data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]

kk <- enrichKEGG(gene         = gene,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
head(kk)

#====非模式物种的KEGG富集分析enricher======

#第一步：准备数据
# ===1-关键基因的列表：pull将一个二维数据转化成一个一维数据，即直接返回一个向量
gene <- filter(de_result, 
               abs(log2FoldChange) > 1 & padj < 0.05) %>% 
  pull(id)
# ===2-TERM与基因对应关系信息
pathway2gene <- select(gene_mapper,
       GID = Gene_Id,
       KO = Ko,
       Pathway = Pathway) %>%
  dplyr::select(Pathway, GID) %>%
  separate_rows(Pathway, sep=",", convert = F) %>%
  filter(str_detect(Pathway, "ko")) %>%
  mutate(Pathway = str_remove(Pathway, "ko"))

#===3-TERM2name用于获得kegg中path对应的name（描述）
library(magrittr)
get_path2name <- function(){
  keggpathid2name.df <- clusterProfiler:::kegg_list("pathway")
  keggpathid2name.df[,1] %<>% gsub("path:map", "", .)
  colnames(keggpathid2name.df) <- c("path_id","path_name")
  return(keggpathid2name.df)
}
pathway2name <-get_path2name()
head(pathway2name)

#===4-所有基因的FC信息
gene_list <- de_result$log2FoldChange
names(gene_list) <- de_result$id #有名字的向量
sort(gene_list, decreasing = T)


#第二步：KEGG富集分析
library(clusterProfiler)
de_ekp <- enricher(gene,
                   TERM2GENE = pathway2gene,
                   TERM2NAME = pathway2name,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05)

de_ekp_df <- as.data.frame(de_ekp)
head(de_ekp_df)
save(de_ekp, de_ekp_df, file = "Data/STEP4_KEGG_enrichment.RData")

library(enrichplot)
barplot()


barplot(de_ekp, showCategory = 15)
dotplot(de_ekp, showCategory = 15)
cnetplot(de_ekp, 
         foldChange = gene_list,
         #如果是个圈的话，这里可以选择节点标签 category | gene | all | none
         node_label = "category",
         circular = T,
         colorEdge = T)
emapplot(de_ekp, showCategory = 10, pie = 'count')

library(pathview)
?pathview





  