library(WGCNA)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(ggthemes)

#WGCNA的核心是什么：
#怎样确定连个基因是否有关系呢？
#相关系数矩阵=> 相关系数幂矩阵=> TOM矩阵(考虑到基因之间的间接的关系)


gene_expr <- read.csv("E:/Rhome/RNA-Seq-pipeline/Data/WGCNA/LiverFemale.csv", row.names=1)
sample_infor <- read.csv("E:/Rhome/RNA-Seq-pipeline/Data/WGCNA/ClinicalTraits.csv", row.names=1)
gene_infor <- read.csv("E:/Rhome/RNA-Seq-pipeline/Data/WGCNA/GeneAnnotation.csv") #不需要变成dataframe

dim(gene_expr)
datExpr0 <- t(gene_expr) #WGCNA要求基因为列

#缺失数据及无波动数据过滤
#一个样本中一半以上的基因都没有表达量去掉，一般芯片这样处理
#普通测序项目中不需要
gsg <- goodSamplesGenes(
  datExpr0,
  minFraction = 1/2 #基因缺失数据比例阈值
)
datExpr <- datExpr0[gsg$goodSamples, gsg$goodGenes]
dim(gene_expr)

#通过聚类，查看是否有明显异常样本，如果有需要剔除
plot(hclust(dist(datExpr)),
     cex.lab = 1.5,
     cex.axis = 1.5,
     cex.main = 2,
     cex = 0.65)
 
#WGCNA作者建议将所有基因进行分析而不是只取差异基因
#2W个基因需要=> 16G ram
#3w个基因需要=> 32G ram
#8k-1w个基因需要 => 4G ram

#样本信息表和样本的表达矩阵的行名一定要对上，最后要对两个
#表格算相关系数
#注意：datTraits中的数据必须全部为数字
#注意：datExpr的数据必须是进行标准化后的(TPM)，wgcna就是用来算相关系数
datExpr[1:4, 1:4]
datTraits <- sample_infor
datTraits[1:4, 1:4]


#第一步：寻找最佳beta值
if(T){
  enableWGCNAThreads(nThreads = 3) #开启多线程

  stf <- pickSoftThreshold(
    datExpr,
    powerVector = 1:20, #尝试1到20
    networkType = "unsigned")
  
  save(stf, file = "Outcomes/WGCNA/stf.RData")
  
  fig_power1 <- ggplot(data = stf$fitIndices, aes(x = Power, y = SFT.R.sq)) +
    geom_point() +
    geom_text_repel(aes(label = Power)) +
    geom_hline(yintercept = 0.85, color = "red") +
    labs(title = "Scale independence",
         x = "Soft Threshold(power)",
         y = "Scale Free Topology Model fit, signed R^2") +
    theme_bw() +
    theme(plot.title = element_text(size = 18, hjust = 0.5))
  
  fig_power2 <- ggplot(data = stf$fitIndices, aes(x = Power, y = mean.k.)) +
    geom_point() +
    geom_text_repel(aes(label = Power)) +
    labs(title = "Mean connectivity", 
         x = "Soft Threshold(power)",
         y = "Mean Connectivity") +
    theme_bw() +
    theme(plot.title = element_text(size = 18, hjust = 0.5))
  
  plot_grid(fig_power1, fig_power2)}

# 第二步：构建网络
if(T){
  net <- blockwiseModules(
    #输入数据
    datExpr,
    #1.计算相关系数，基因与基因之间的相关性
    corType = "pearson",
    #2.计算邻接矩阵
    power = stf$powerEstimate,
    networkType = "unsigned",
    #3.计算TOM矩阵
    TOMType = "unsigned",
    saveTOMs = T,
    saveTOMFileBase = "blockwiseTOM",
    #4.聚类并划分模块
    deepSplit = 2, #0|1|2|3|4,值越大模块就越多越小
    minModuleSize = 30, # 保留的模块的最小基因数
    #5.合并相似模块
    #计算模块特征向量(module eigengenes,MEs),
    #计算MEs与datTrait之间的相关性
    #对距离小于mergeCutHeight(1-cor)的模块进行合并
    mergeCutHeight = 0.25,#相关性大于0.75的合并 
    
    numericLabels = F,#数字命名模块
    nThreads = 0, #使用所有可用线程
    maxBlockSize = 100000) #输入的最多基因数
  table(net$colors)}

if(T){
  library(tidyverse)
  wgcna_result <-data.frame(gene_id = names(net$colors),
                            module = net$colors) %>%
    left_join(gene_infor, by = c("gene_id" = "gene_id"))
  head(wgcna_result)
  save(wgcna_result, file = "Outcomes/WGCNA/wgcna_result.RData")
}
 
# 模块的可视化
plotDendroAndColors(
  dendro = net$dendrograms[[1]],
  colors = net$colors,
  groupLabels = "Module colors",
  dendroLabels = F,
  addGuide = T)

#模块特征向量(module eigengenes,MEs)，简单说我们需要计算出：每个模块和性状的想关系数，但是模块中
#中的基因很多，相关性的计算只能两两计算，因此wgcna使用pac，得出模块的pc1，然后求pac1与性状之间
#的相关系数，模块特征向量就是pc1

#第三步：计算模块（模块特征向量）与性状之间的相关性
moduleTraitCor <- cor(
  net$MEs, datTraits,
  use = "p",
  method = "pearson" #注意相关系数的计算方式
)

#计算pvalue
moduleTraitPvalue <- corPvalueStudent(
  moduleTraitCor,
  nrow(datExpr))

if(T){#可视化模块与性状的相关性
  sizeGrWindow(10,6)
  # 连接相关性和 pvalue
  textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                      signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) <- dim(moduleTraitCor)
  # heatmap 画图
  par(mar = c(6, 8.5, 3, 3))
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = names(datTraits),
                 yLabels = names(net$MEs),
                 ySymbols = names(net$MEs),
                 colorLabels = FALSE,
                 colors = greenWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))}

if(T){#另外两种可视化模块与性状相关性的方法,如果性状较少建议使用
  MET <- orderMEs(cbind(net$MEs, dplyr::select(datTraits, weight_g)))
  plotEigengeneNetworks(
    multiME = MET, 
    setLabels = "Eigengene dendrogram", 
    plotDendrograms = TRUE, 
    plotHeatmaps = FALSE,
    colorLabels = TRUE,
    marHeatmap = c(3,4,2,2))}
if(T){#另外两种可视化模块与性状相关性的方法,如果性状较少建议使用
  plotEigengeneNetworks(
    multiME = MET, 
    setLabels = "Eigengene dendrogram", 
    plotDendrograms = FALSE, 
    plotHeatmaps = TRUE,
    colorLabels = TRUE,
    marHeatmap = c(8,8,2,2))}

# 提取模块中的基因
getgene <- function(module = "color", threshold = 0.2, power = 6, threshold1 = 0.2){
  my_modules <- module
  m_wgcna_result <- filter(wgcna_result, module %in% my_modules)
  m_datExpr <- datExpr[, m_wgcna_result$gene_id]
  m_TOM <- TOMsimilarityFromExpr(
    m_datExpr,
    power = power,
    networkType = "unsigned",
    TOMType = "unsigned")
  dimnames(m_TOM) <- list(colnames(m_datExpr), colnames(m_datExpr))
  cyt <- exportNetworkToCytoscape(
    m_TOM,
    edgeFile = "Outcomes/WGCNA/CytoscapeInput-network.txt",
    weighted = TRUE,
    threshold = threshold1)}

getgene(module = "blue") #提取蓝色模块中的基因

#找hub基因，用网络连接度这种方式来找，谁的degree高谁重要，
# 网络degree
# 1.连接数
# 2.中介连接度



