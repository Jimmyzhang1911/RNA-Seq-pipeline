load("Data/STEP2_de_result.RData")
load("Data/STEP1_input.RData")
source("utils.R")
library(tidyverse)
library(pheatmap)



#火山图：看差异基因在所有基因的大概分布，看不出每个差异基因在样本中的具体表达量，热图可以看出来
top_de_expr <- dplyr::slice(de_result, c(1:20)) %>% 
  select(-c(2:7)) %>%
  column_to_rownames(var = "id")
write.csv(top_de_expr, file = "Outcomes/top20_de_gene.csv", quote = F)

pheatmap(top_de_expr) #问题原因：差异太大

#解决方法一：用log拉近距离，+1是因为有的数据是0。
#只是数据的波动范围缩小，数据原有的大小差异依然存在
#例如也就是说同一样品中，基因的表达差异依然可以看出来
colr <- list(
  stage = c(S1 = "#7FFF00", S2 = "#006400", S3 = "#FFD700", S4 = "#FF0000")
)
pdf("Figures/heatmap.pdf")
plot = pheatmap(log10(top_de_expr + 1),
         #不对行聚类
         #cluster_rows = F,
         #show_colnames = F,
         annotation_col = select(sample_info, 2),
         annotation_colors = colr,
         #color = colorRampPalette(c("green", "white", "red"))(100),
         cutree_rows = 3,
         cutree_cols = 2,
         angle_col = 45) 
dev.off()
#解决方法二：标准化(两层含义)：scale两件事情同时做
# 标准化：把数据变化的幅度，归一化到一个统一范围内
# 中心化：把数据的平均值挪到0
#不同基因在不同样本不具可比性
pheatmap(top_de_expr, 
         scale = "row") #按行标准化
colorRampPalette(c("blue","red"))[200]