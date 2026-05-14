library(dplyr)
library(dbplyr)
library(stringr)
library(ggplot2)
library(viridis)
library(clusterProfiler)
library(ChIPseeker)
library(GenomicFeatures)
library(GenomicRanges)
library(chromVAR) ## For FRiP analysis and differential analysis
library(DESeq2) ## For differential analysis section
library(ggpubr) ## For customizing figures
library(corrplot) ## For correlation plot
library(biomaRt)
library(curl)
library(org.Dm.eg.db)

getwd()
#setwd("L:\\E75\\CUT&Tag part")
setwd("I:\\PPP1CA&PPP1CC\\02BulkATACseq\\10Chipseeker")

txdb <- makeTxDbFromGFF("Drosophila_melanogaster.BDGP6.32.109.gtf",
                        format="gtf")    #可以使用gtf和gff3
						
						
keytypes(txdb)    #感兴趣的话，可以用以下方法探索txdb都包含了什么内容
keys(txdb)

#读入单个summits文件
peaks <- readPeakFile("macs2_E75_Flag_peak_q0.05_summits.bed")
peaks <- readPeakFile("macs2_NY_N_peak_q0.05_summits.bed")
peaks <- readPeakFile("macs2_NY_Y_peak_q0.05_summits.bed")
peaks <- readPeakFile("macs2_ENY_N_peak_q0.05_summits.bed")
peaks <- readPeakFile("macs2_ENY_Y_peak_q0.05_summits.bed")
#结构注释
peakAnno <- annotatePeak(peaks, TxDb=txdb, tssRegion=c(-2000, 2000))

#最后将我们的注释结果转为数据框，便于查看
df <- as.data.frame(peakAnno)
head(df)
#将注释到的基因提取出来（第14列），用于后续功能分析
gene <- df[,14]
###对基因进行注释-获取gene_symbol
gene.df <- bitr(gene, fromType = "ENSEMBL", 
                toType = c("SYMBOL","ENTREZID"),
                OrgDb = org.Dm.eg.db) 
# geneID    需要转换的ID
# fromType  当前ID类型
# toType    转换成什么ID，使用keytypes()查看有哪些类型
# OrgDb     注释数据库
# 注意相同的FBgn号会注释出不同的基因，所以以下代码会报错
#row.names(gene.df) <- gene.df$ENSEMBL重复 FBgn0004828对应 His3.3B和His3.3A
#row.names(gene.df) <- gene.df$SYMBOL重复 His2B:CG33872对应 FBgn0053910和FBgn0053872
write.csv(gene.df,file="peakAnno_NY_Y.csv",sep="/t")


a <- read.table(file = "Venn_ENY_NY_N&Y_712_962交集情况.txt",sep = "\t",header = T)
gene.df2 <- bitr(a[,1], fromType = "ENSEMBL", 
                toType = c("SYMBOL","ENTREZID"),
                OrgDb = org.Dm.eg.db) 
write.csv(gene.df2,file="Anno_ENY_NY_N&Y_79.csv",sep="/t")

a <- read.table(file = "Venn_Only ENY_N.txt",sep = "\t",header = T)
gene.df2 <- bitr(a[,1], fromType = "ENSEMBL", 
                toType = c("SYMBOL","ENTREZID"),
                OrgDb = org.Dm.eg.db) 
write.csv(gene.df2,file="Anno_Venn_Only ENY_N_712.csv",sep="/t")

a <- read.table(file = "Venn_Only ENY_Y.txt",sep = "\t",header = T)
gene.df2 <- bitr(a[,1], fromType = "ENSEMBL", 
                toType = c("SYMBOL","ENTREZID"),
                OrgDb = org.Dm.eg.db) 
write.csv(gene.df2,file="Anno_Venn_Only ENY_Y_962.csv",sep="/t")


#一次也可以读入多个summits文件，使用list存储，然后使用lapply注释
files = list(E75_Flag = ("macs2_E75_Flag_peak_q0.05_summits.bed"), 
				 NY_N = ("macs2_NY_N_peak_q0.05_summits.bed"), 
				 NY_Y = ("macs2_NY_Y_peak_q0.05_summits.bed"),
				 ENY_N = ("macs2_ENY_N_peak_q0.05_summits.bed"), 
				 ENY_Y = ("macs2_ENY_Y_peak_q0.05_summits.bed"))
peakAnnoList <- lapply(files, 
                       annotatePeak,
                       TxDb=txdb,
                       tssRegion=c(-2000, 2000))
plotAnnoBar(peakAnnoList)
plotDistToTSS(peakAnnoList)


#注释完，进行可视化，多种图可供选择
plotAnnoBar(peakAnno)
plotDistToTSS(peakAnno)
vennpie(peakAnno)
plotAnnoPie(peakAnno)
#install.packages("ggupset")
library(ggupset)
upsetplot(peakAnno)
#install.packages("ggimage")
library(ggimage)
upsetplot(peakAnno, vennpie=TRUE)




