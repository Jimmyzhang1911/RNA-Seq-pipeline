load("Data/STEP1_input.RData")

library(tidyverse)
#导入差异分析结果，这里不希望有行名，因为后面要用tidyvers处理 
#select(-A, -B) #反选
de_result <- read.table("Data/genes.counts.matrix.KID_S1_vs_KID_S4.DESeq2.DE_results",
                 header = T)

#过滤列select, 过滤行filter, mutatet新添加一列，left_join注意用=连接共同的列名
#得到差异分析的总表格涵基因表达量，上下调，pvalue，log2FoldChange，基因信息等等
#这个表保存了很多数据很关键
de_result <- mutate(de_result, direction = ifelse(
  padj > 0.05, "ns", ifelse(
    abs(log2FoldChange) < 1, "ns", ifelse(
      log2FoldChange >= 1, "up", "down")))) %>%
  left_join(gene_infor, by = c("id" = "Gene_Id")) %>%
  left_join(rownames_to_column(gene_expr, var = "id"), by = "id") %>%
  select(-c(2:4,6:7)) %>%
  arrange(desc(abs(log2FoldChange)))

#分组统计,n()数一下每组有多少行
table(de_result$direction)
group_by(de_result, direction) %>% 
  summarise(count = n())

save(de_result,file="Data/STEP2_de_result.RData")

load("Data/STEP2_de_result.RData")
#差异分析两张图要放，一个是火山图（找一些基因知道他们的差异情况）
#一个是heatmap（每个基因表达量的变化）
#火山图：
library(EnhancedVolcano)
EnhancedVolcano(de_result,
                #一般使用Gene_Symbol,但是非模式生物的不是很丰富所以先用id
                lab = de_result$id, 
                x = "log2FoldChange",
                y = "padj",
                FCcutoff = 2,
                pCutoff = 0.01)


#ggplot2绘制火山图
library(ggplot2)
library(ggrepel)

#数据，坐标
de_result <- mutate(de_result,
       direction = factor(direction, levels = c("down", "ns", "up")))

top_de <- filter(de_result, abs(log2FoldChange) > 2 & -log10(padj) > 100)
top_de2 <- filter(de_result,
                 id == c("HF01768", "HF12718")) # id %in% c("HF01768", "HF12718")

my_palette <-c("#00008B", "#708090", "#8B0000")


p <- ggplot(data = de_result, aes(x = log2FoldChange, y = -log10(padj))) + 
  geom_point(aes(color = direction, 
                 size = abs(log2FoldChange))) +
  geom_label_repel(data = top_de, aes(label = id)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(1, -1), linetype = "dashed") +
  scale_size(range = c(0.1, 1.5)) +
  scale_color_manual(values = my_palette,
                     labels = c("up 1221 genes", "ns 20925 genes", "down 1080 genes"),
                    name = "S1 vs S2") +
  labs(x = "log2 FOldChange",
       y = "-log10(pvalue)",
       title = "DE Volcano plot",
       size = "log2FOldChange") +
  ylim(c(0,200)) +
  guides(size = F) +
  theme_bw() +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.position = c(0.14, 0.84))
p

#分步解释
p <- ggplot(data = de_result, aes(x = log2FoldChange, y = -log10(padj))) + 
  # 图形为点,以及点的属性，ase=>映射，alpha透明度，shape性状不应该添加到映射中
  geom_point(aes(color = direction, 
                 size = abs(log2FoldChange))) +
  # 添加标签,size设置大小
  geom_label_repel(data = top_de, aes(label = id)) +
  # 添加线条,hline水平，size线的宽度，color线的颜色
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  # 添加垂直测线
  geom_vline(xintercept = c(1, -1), linetype = "dashed") +
  # 设置属性波动的范围
  scale_size(range = c(0.1, 1.5)) +
  # 设置颜色，可将数据改为factor，然后对应，标度
  scale_color_manual(values = my_palette,
                     labels = c("up 1221 genes", "ns 20925 genes", "down 1080 genes"),
                     name = "S1 vs S2") +
  # 设置坐标轴名称
  labs(x = "log2 FOldChange",
       y = "-log10(pvalue)",
       title = "DE Volcano plot",
       size = "log2FOldChange") +
  # 设置坐标轴范围
  ylim(c(0,200)) +
  # 修改图例
  guides(size = F) +
  # 设置主题
  theme_bw() +
  # 修改主题,水平上看左对齐右对齐hjust，垂直上用vjust，（0，0.5，1）=>左中右
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        legend.background = element_blank(),
        legend.key = element_blank(), #去掉圆点后面的颜色
        legend.position = c(0.14, 0.84))
ggsave(p,filename = "Figures/DE_volcano_plot.pdf",width = 12,height = 9)    

