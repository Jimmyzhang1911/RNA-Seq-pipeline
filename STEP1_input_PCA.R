source("utils.R")
library(tidyverse)
library(readr)
library(pheatmap)
library(PCAtools)

#===导入数据，数据结构：有行名的叫dataframe，没有行名的叫table
#基因表达矩阵为TMM标准化后的可做PCA
gene_expr <- read.table("Data/genes.TMM.EXPR.matrix", header = T, row.names = 1) 
#样本信息表
sample_info <- read.table("Data/sample_info.txt", header = T, row.names = 1) 
#基因信息表，数据结构用table，select选择列，可同时换列的名字
gene_mapper <- read_delim("Data/query_seqs.fa.emapper.annotations", 
                         "\t", escape_double = FALSE, col_names = FALSE, 
                         comment = "#", trim_ws = TRUE) %>%
  select(Gene_Id = X1, 
         Gene_Symbol = X6,
         Go = X7,
         Ko = X9,
         Pathway = X10,
         COG = X21,
         Gene_Name = X22)
rownames(sample_info) <- sub("LD", "Rep", rownames(sample_info)) #名称不对称，替换一下，后续PCA要用
gene_infor <- dplyr::select(gene_mapper, Gene_Id, Gene_Symbol, Gene_Name) %>%
  filter(!is.na(Gene_Name))

save(gene_expr, sample_info, gene_mapper, gene_infor, file="Data/STEP1_input.RData")
load("Data/STEP1_input.RData")
load("Data/STEP2_de_result.RData")

#===STEP1-先做样本相关性分析(再做差异分析)
#算法：pearson | kendall | spearman 再method中选择
# pearson：线性相关，比如两个基因的相关性，样本之间相关系数
# kendall：离散型相关，分类型，比如性别相关
# spearman：等级相关，比如那些基因与肿瘤的分期相关，肿瘤的分期，一期二期三期四期，他们之间是等级关系不是线性关系
sample_cor <- round(cor(gene_expr, method = "pearson"), digits = 2) #算出样本间的相关系数，保留两位小数 
pheatmap(sample_cor) #画个热图简单看一下，行列皆为样本描述样本之间的相关性

#===STEP2-样本聚类
#计算距离矩阵，任意两两样本之间的离的有多远，相似度有多高，dist函数求的是行的距离
#欧几里得eudidean，铁米雪夫maxmum，一般用欧几里得不用改
sample_dist <- dist(t(gene_expr))
#第二步：聚类
#算法：层次聚类法，k-means划分聚类，一般用层次聚类法
sample_hc <- hclust(sample_dist)
plot(sample_hc)  
   
#===STEP3-PCA主成分分析法
#目的：之前每一个基因都描述了样本之间的差异，PCA将几万个基因=>几个对差异贡献最大的几个成分
p <- pca(gene_expr, metadata = sample_info, removeVar = 0.1) #bioinfor中将样本信息表称为metadata
p$components #看一眼降维成几个主成分
p$rotated #每个样本与每个主成分之间的相关性
p$loadings #每个基因对每个主成分的贡献度

#每个主成分都描述了样本之间的差异，有的主成分能拉开样本之间的差异，
#有的拉不开，PC1对样本差异的解释度最高
screeplot(p)

#每个点是一个样本，可以看出PC1刚好把时期区分开来，
#如果想研究和时期相关的基因就可以去找对PC1贡献最大的基因
biplot(p, x = "PC1", y = "PC2", 
       colby = "strain",
       colkey = c("KID" = "red", "BLO" = "blue"),
       shape = "stage",
       legendPosition = "right") 
 
# 可以看出品种之间的差异，可以多看几个主成分
biplot(p, x = "PC1", y = "PC3", 
       colby = "strain",
       colkey = c("KID" = "red", "BLO" = "blue"),
       shape = "stage",
       legendPosition = "right") 

plotloadings(p)#那些基因对相应的主成分贡献度大


#一般项目中，组类差异小，组间差异大
#批次效应：用一批的样本聚到一起
#看一眼批次效应：中位数是不是接近，TMM可以消除部分批次效应，
#去除批次效应主要是将各个样中基因表达的中位数拉平，可用SVA这个R包，
boxplot(log10(gene_expr + 1)) 

