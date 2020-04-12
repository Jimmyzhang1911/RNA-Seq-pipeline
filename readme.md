[TOC]
# 转录组流程：
## 软件安装：
## 数据准备：
测序数据：sample.fastq
基因组序列：genome.fasta
基因注释文件：gene.gtf
蛋白序列：proteins.fasta
样本信息表：samples.txt
## STEP1_比对：
**hisat2**: sample.fastq + genome.fasta => sample.bam (每个序列比对到基因组的位置)
step1: 构建index:    `hisat-build`   sample.fastq => index/genome
step2: 比对`hisat2`sample.fastq + index => sample.sam
step3: 压缩和排序`samtools sort`sample.sam => sample.bam 
step4: 构建bam文件索引`samtools index`sample.bam => sample.bam.bai(IGV对比对结果进行可视化需要)
```
#使用awk和samples.txt批量生成脚本
$ hisat-build ../ref/genome.fasta ../ref/genome 1>hisat2-build.log 2>&1 #构建index
$ hisat2 --new-summary -p 2 -x ../ref/ht_index/genome -U /home/ug0007/jimmy/rna-seq-pipeline/clean_data/BLO_S1_Rep1.fastq.gz -S ../bam/BLO_S1_Rep1.sam --rna-strandness R 1> ../bam/BLO_S1_Rep1.log 2>&1 #比对
$ samtools sort -o ../bam/BLO-S1_Rep1.bam ../bam/BLO_S1_Rep1.sam & #压缩和排序
$ samtools index ../bam/BLO-S1_Rep1.bam & #构建bam的文件索引
$ samtools flagstat ../bam/BLO-S1_Rep1.bam # 查看samtools统计的比对率
```

>`samtools view sample.bam | less -S`查看bam文件 
`--new-summary` :hisat2的前身tophat会产生日志文件，如果想显示新的日志文件格式就输入

###IGV简单使用
genome.fasta和genes.gtf在IGV中创建新的genome

因为用的是测试数据所以效果不是特别好，一般都比对到外显子上，每一个灰色的条就是一条reads，方向都往一个方向，因为用`hisat2`比对的时候设置了参数`--rna-strandness`：链特异性文库则按同一个方向比对。reads上彩色位点是与基因组上不一样的位点，如果只有某一条不一样则大概率可能是测序错误也有可能是RNA-Edit，如果一半一样一半不一样则大概率这个位点为一个杂合。DNA突变，DNA比对某个位点大部分与参考基因组不一样。**个体重测序的本质是基因分型，群体重测序的本质是等位基因的频率**




## STEP2_定量（每个基因的表达量）:
**subread** 中的`featurecounts`:   sample.bam（比对到的位置） + gene.gtf（每个基因的位置） => sample.count
> sunbread 有R和linux两个版本：
- Linux:`conda install bioconductor-rsubread`
- R:`biocManager::install("subread")`

> R的软件包来源：
- cran `insatll.packages("ggplot2")` `conda install r-ggplot2`
- bioconductor `biocManager::install("subread")` `conda install bioconductor -rsubread`

> 有参的定量软件分为两类：基于比对的（featurecounts,htseq），不基于比对的而是基于K-mer频率（salmon ）
转录本的拼接组装：依据测试数据比对到基因组上的结果和已有基因结构的信息对已有的基因结构进行修正，优化和延长不用stringtie的原因：stringtie实现两种功能：转录本的拼接和计算表达量
```
$ Rscript run-featurecounts.R -b ./bam/KID_S1_Rep1.bam -g ./ref/genes.gtf -o ./counts/KID_S1_Rep1 1> ./counts/KID_S1_Rep1.log 2>&1 #定量，计数
$ ls
KID_S1_Rep1.count
KID_S1_Rep1.log
$ cat KID_S1_Rep1.log
Assigned        298072 #比对到基因组的reads数     
Unassigned_Unmapped     130267 #没有比对上到基因组的reads数
Unassigned_Read_Type    0
Unassigned_Singleton    0
Unassigned_MappingQuality       0
Unassigned_Chimera      0
Unassigned_FragmentLength       0
Unassigned_Duplicate    0
Unassigned_MultiMapping 0
Unassigned_Secondary    0
Unassigned_NonSplit     0
Unassigned_NoFeatures   89007 #比对到基因组上但是没有基因结构，如果太大说明基因组注释不到位
Unassigned_Overlapping_Length   0
Unassigned_Ambiguity    554
$ head -5 KID_S1_Rep1.count
gene_id counts  fpkm    tpm
HF32875 0       0       0
HF32876 15      26.0878236356267        24.6606288040338
HF32877 3       2.08508024831671        1.97101110267169
HF32878 25      30.5212347119868        28.8514998562477
```
> 比对上但是没有基因结构的地方算不算新基因？1.是一个新基因，这个外显子来自另一个基因；2.新的转录本，只是这个基因的一部分；3.基因结构优化，之前gtf里面的基因结构有问题
## STEP3_合并成矩阵
`abundance_estimates_to_matrix.pl`
**sample.count**  =>  reads count矩阵 (**gene_counts.matrix**) + 标准化后的矩阵(**TPM.matrix**)
```
$ vim merge_result.sh
#!/bin/bash
ls ./counts/*.count > ./matrix/genes.quant_files.txt
perl abundance_estimates_to_matrix.pl --est_method featureCounts --quant_files ./matrix/genes.quant_files.txt --out_prefix ../matrixgenes

bash merge_result.sh
$ ls
genes.counts.matrix #标准化之前的counts用来做差异分析，因为DESeq2和edgeR自带标准化
genes.TMM.EXPR.matrix #经过TMM（样本间的标准化）标准化过后（消除部分批次xiaoying），用来做共表达分析，PCA分析，heatmap等等
genes.TPM.not_cross_norm #经过TPM（样本内的标准化）标准化过后，消除基因长度的影响，不同样本测序深度的影响
genes.TPM.not_cross_norm.TMM_info.txt
```
## STEP4_差异表达分析
> 组间差异远大于组内差异则认为是差异基因

 reads count矩阵 + 样本信息表(samples.txt) + 比对设计表(contrasts.txt)
用到`Trinity`中的`run_DE_analysis.pl`脚本
```
$ conda install -n rna-seq trinity
$ which -a run_DE_analysis.pl
/home/ug0007/miniconda3/envs/rna-seq/bin/run_DE_analysis.pl #当前环境变量的
/pub/anaconda3/bin/run_DE_analysis.pl #公共环境变量中的
#当前的后来用不了因为加载不了edgeR，尝试在当前环境用conda装网不好，没能成功于是用的公共环境中的

$ vim run_DE_analysis.run
#!/bin/bash
/pub/anaconda3/bin/run_DE_analysis.pl \    #公共环境变量中的
        --matrix ../matrix/genes.counts.matrix \ 
        --method DESeq2 \
        --samples_file ../clean_data/samples.txt \
        --contrasts contrasts.txt #比对设计表
$ ls 
genes.counts.matrix.KID_S1_vs_BLO_S1.DESeq2.DE_results  #结果 |LFC(差异的倍数)| > 1 && padj < 0.05
genes.counts.matrix.KID_S1_vs_BLO_S1.DESeq2.Rscript #脚本，本质是用perl生成R脚本然后再linux上运行
genes.counts.matrix.KID_S1_vs_BLO_S1.DESeq2.DE_results.MA_n_Volcano.pdf #图

#找出差异基因
sed '1d' genes.counts.matrix.KID_S1_vs_BLO_S1.DESeq2.DE_results | less -S #删除表头
sed '1, 10d' genes.counts.matrix.KID_S1_vs_BLO_S1.DESeq2.DE_results | less -S #删除第一行到第十行
sed '1d' genes.counts.matrix.KID_S1_vs_BLO_S1.DESeq2.DE_results | awk 'sqrt($5*$5) > 1 && $9 < 0.05 {print $1}' | wc
# 打印重要的几列，按第二列排序，-k 2按字母表，-k 2n按数字大小
sed '1d' genes.counts.matrix.KID_S1_vs_BLO_S1.DESeq2.DE_results | awk 'sqrt($5*$5) > 1 && $9 < 0.05 {print $1"\t"$5"\t"$9}' | sort -k 2n
```
## STEP5_功能注释：eggNOG-mapper的使用
eggNOG是一个基因家族数据库，数据量没有Nr，interpro的多但是一般够用
网站：http://eggnog5.embl.de/#/app/home
提交CDS或者蛋白的FASTA文件注意查收邮件，需要通过邮件中的链接启动项目：




