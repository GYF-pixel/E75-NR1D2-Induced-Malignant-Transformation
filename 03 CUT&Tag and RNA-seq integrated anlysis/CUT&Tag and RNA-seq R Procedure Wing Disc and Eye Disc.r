library(devtools)
library(biomaRt)
library(curl)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(pheatmap)
library(RColorBrewer)
library(tidyr)
library(stringr)
library(edgeR)
library(ggsci)
library(cowplot)
library(tidyverse)
library(ggunchull)
library(SCENIC)
library(scales)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
#library(xlsx)
setwd("F:\\E75\\CUT&Tag part\\07RNAseq Combined")                    #设置工作目录 EyeDisc

##Step 2 DEGs identification by edgeR
library(devtools)
library(biomaRt)
library(curl)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(pheatmap)
library(RColorBrewer)
library(tidyr)
library(stringr)
library(edgeR)
library(ggsci)
library(cowplot)
library(tidyverse)
library(ggunchull)
library(SCENIC)
library(scales)
#library(xlsx)
library(dbplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
#setwd("L:\\E75\\RNA_seq part\\04Expression\\WingDisc")                    #设置工作目录 EyeDisc 替换为 WingDisc
setwd("L:\\E75\\RNA_seq part\\04Expression\\EyeDisc")     
setwd("F:\\E75\\CUT&Tag part\\07RNAseq Combined")

##foldChange=0.5849625    foldChange=log2(FC=1.5)=0.5849625
foldChange=0.5849625
padj=0.05
#foldChange=2
padj=0.05

##手动删除all_gene_count.tsv没有数据对应的行 
rt=read.csv("newdata_filter.csv")              #改成自己的文件名
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
data = exp
data.filter = data[rowMeans(data) > 1 & rowSums(data) > 15, ]   #数据过滤

# condition table
colnames(data.filter)
group_info = c(rep("ENY",3), rep("NY",3), rep("WT",3))

#DGEList constructor
dge.list.obj <- DGEList(counts = data.filter, group = group_info)
dge.list.obj
# Normalization method: "TMM","TMMwsp","RLE","upperquartile","none"
dge.list.obj <- calcNormFactors(dge.list.obj,method = "RLE") # DESeq2, cuffdiff
dge.list.obj$samples
#plotMDS(dge.list.obj)
# make design matrix
design.mat <- model.matrix(~group_info)
# estimate dispersion
dge.list.obj <- estimateDisp(dge.list.obj,design.mat)
dge.list.obj$common.dispersion
dge.list.obj$tagwise.dispersion
	### Equal to the following steps
	#### 1st common dispersion
	#dge.list.obj <- estimateCommonDisp(dge.list.obj)
	#### 2nd tagwise dispersion
	#dge.list.obj <- estimateTagwiseDisp(dge.list.obj)
	#### plot dispersion
	#plotBCV(dge.list.obj, cex = 0.8)
	#### plot var and mean
	#plotMeanVar(dge.list.obj, show.raw=TRUE, show.tagwise=TRUE, show.binned=TRUE)

# test with likelihood ratio test
#1.  GLM general linear model
	#fit <- glmFit(dge.list.obj, design.mat)
	#lrt <- glmLRT(fit, coef=2)
	#DEGs.res.lrt <- as.data.frame(topTags(lrt,n=nrow(count_df.filter),sort.by = "logFC"))
#2.  exactTest
#dge.list.res <- exactTest(dge.list.obj, pair = c("WT","NY"))
#dge.list.res <- exactTest(dge.list.obj, pair = c("WT","ENY"))
dge.list.res <- exactTest(dge.list.obj, pair = c("NY","ENY"))
#topTags 按不同的标准来排列差异数据
topTags(dge.list.res)
ordered_tags <- topTags(dge.list.res, n=1000000, sort.by = "logFC")
head(ordered_tags)

DEGs.res <- as.data.frame(topTags(dge.list.res,n=nrow(data.filter),sort.by = "logFC"))

# 输出相关的文件
allDiff=DEGs.res[is.na(DEGs.res$FDR)==FALSE,]
diff=allDiff
###对基因进行注释-获取gene_symbol
library("biomaRt")
library("curl")
mart <- useDataset("dmelanogaster_gene_ensembl", useMart("ensembl"))
my_ensembl_gene_id<-row.names(diff)
study_symbols<- getBM(attributes=c('ensembl_gene_id','external_gene_name',"description"), filters = 'ensembl_gene_id', values = my_ensembl_gene_id, mart = mart)
head(study_symbols)
ensembl_gene_id<-rownames(diff)

diff <- cbind(ensembl_gene_id,diff)
colnames(diff)[1]<-c("ensembl_gene_id")
#全部基因的edger输出结果--diff_name
diff_name <- merge(diff,study_symbols,by="ensembl_gene_id")
write.table(diff_name,file="Part1_edgerOut.xls",sep="\t",quote=F)
#所有的差异基因--diff_nameSig
diff_nameSig = diff_name[(diff_name$FDR < padj & (diff_name$logFC>foldChange | diff_name$logFC<(-foldChange))),]
write.table(diff_nameSig, file="Part1_diff_nameSig.xls",sep="\t",quote=F)
nrow(diff_nameSig)
#所有的上调基因--diff_nameUp
diff_nameUp = diff_name[(diff_name$FDR < padj & (diff_name$logFC>foldChange)),]
write.table(diff_nameUp, file="Part1_diff_nameup.xls",sep="\t",quote=F)
nrow(diff_nameUp)
#所有的下调基因--diff_nameDown
diff_nameDown = diff_name[(diff_name$FDR < padj & (diff_name$logFC<(-foldChange))),]
write.table(diff_nameDown, file="Part1_diff_namedown.xls",sep="\t",quote=F)
nrow(diff_nameDown)
#所有基因normalized counts table--normalizeExp
normalizeExp=cpm(dge.list.obj)
write.table(normalizeExp,file="Part1_normalizeExp.xls",sep="\t",quote=F,col.names=F)   #输出所有基因校正后的表达值（normalizeExp.txt）

#DEGs normalized counts table--diffExp
diffExp=normalizeExp[diff_nameSig$ensembl_gene_id,]
write.table(diffExp,file="Part1_diffmRNAExp.xls",sep="\t",quote=F,col.names=F)         #输出差异基因校正后的表达值（diffmRNAExp.txt）
nrow(diffExp)

# MA plot
select.sign.gene_id = rownames(dge.list.res)[as.logical(rownames(diffExp))]
#限制y轴要用 ylim=c(-20,20)
MAplot <- plotSmear(dge.list.res, de.tags = rownames(diffExp), cex = 0.7,xlab="Average Log2CPM", ylab="Log2FC",ylim=c(-20,20))

pdf(file="./Part2_MAplot.pdf",width = 6,height = 6)
MAplot <- plotSmear(dge.list.res, de.tags = rownames(diffExp), cex = 0.7,xlab="Average Log2CPM", ylab="Log2FC",ylim=c(-20,20))
dev.off()

#abline(h = c(-0.5849625, 0.5849625), col = "blue")  #添加分割线

# Step 2 Volcano Plot
#画火山图
library("ggplot2")
library("ggpubr")
library("ggthemes")
library("pheatmap")
#在画图之前需要将padj转换成-1*log10，这样可以拉开表达基因之间的差距    -log10(FDR=0.05)=1.30103 log2(FC=1.5)=0.5849625 -log10(FDR=0.01)=2
diff_name2 <- diff_name
rownames(diff_name2)=diff_name2[,1]
diff_name2$log10diff_namepadj <- -log10(diff_name2$FDR)
diff_name2$group <- "nonsignificance"
diff_name2$group[(diff_name2$log10diff_namepadj > 1.30103) & (diff_name2$logFC > foldChange)]="up-regulated"
diff_name2$group[(diff_name2$log10diff_namepadj > 1.30103) & (diff_name2$logFC < -foldChange)]="down-regulated"
table(diff_name2$group)
#Label-FDR最高的20个基因(ensembl_gene_id和external_gene_name可以互换)
diff_name2$label=""
diff_name2 <- diff_name2[order(diff_name2$FDR),]
diff_name2_upgenes_FDR <- head(diff_name2$external_gene_name[which(diff_name2$group=="up-regulated")],20)
diff_name2_downgenes_FDR <- head(diff_name2$external_gene_name[which(diff_name2$group=="down-regulated")],20)
diff_name2_top20genes_FDR <- c(as.character(diff_name2_upgenes_FDR),as.character(diff_name2_downgenes_FDR))
diff_name2$label[match(diff_name2_top20genes_FDR,diff_name2$external_gene_name)] <- diff_name2_top20genes_FDR
#Label-logFC最高的10个基因
diff_name2$logFC_abs=abs(diff_name2$logFC)
diff_name2 <- diff_name2[order(diff_name2$logFC_abs),]
diff_name2_upgenes_logFC <- tail(diff_name2$external_gene_name[which(diff_name2$group=="up-regulated")],10)
diff_name2_downgenes_logFC <- head(diff_name2$external_gene_name[which(diff_name2$group=="down-regulated")],10)
diff_name2_top10genes_logFC <- c(as.character(diff_name2_upgenes_logFC),as.character(diff_name2_downgenes_logFC))
diff_name2$label[match(diff_name2_top10genes_logFC, diff_name2$external_gene_name)] <- diff_name2_top10genes_logFC
write.table(diff_name2,file="Part1_diff_name2.xls",sep="\t",quote=F,col.names=T)   
#画分界线（1.30是-log10(FDR=0.05)的值） log10(FDR=0.05)=-1.30103 log2(FC=1.5)=0.5849625
volcano <- ggscatter(data = diff_name2,x = "logFC",y = "log10diff_namepadj",
	color = "group", 
	palette = c("#2f5688","#BBBBBB","#CC0000"),
	size = 1,
	font.label = c(8, "plain"),
	label = diff_name2$label,
	repel = T,
	xlab="Log2FoldChange",
	ylab="-Log10(FDR)",)+theme_base()+
	geom_hline(yintercept = 1.3,linetype="dashed")+ 
	geom_vline(xintercept = c(-foldChange,foldChange),linetype="dashed")
ggsave("Part2_Volcanoplot.pdf", plot = volcano, width = 11, height = 8) 


# Step 2 Volcano Plot New
#From edgeR results to count_matrix
DEGs.res <- as.data.frame(topTags(dge.list.res,n=nrow(data.filter),sort.by = "logFC"))
# 输出相关的文件
allDiff=DEGs.res[is.na(DEGs.res$FDR)==FALSE,]
diff=allDiff
###对基因进行注释-获取gene_symbol
library("biomaRt")
library("curl")
mart <- useDataset("dmelanogaster_gene_ensembl", useMart("ensembl"))
my_ensembl_gene_id<-row.names(diff)
study_symbols<- getBM(attributes=c('ensembl_gene_id','external_gene_name',"description"), filters = 'ensembl_gene_id', values = my_ensembl_gene_id, mart = mart)
head(study_symbols)
ensembl_gene_id<-rownames(diff)
diff <- cbind(ensembl_gene_id,diff)
colnames(diff)[1]<-c("ensembl_gene_id")
#全部基因的edger输出结果--diff_name
diff_name <- merge(diff,study_symbols,by="ensembl_gene_id")

#准备需要的文件diff_name2--input
diff_name2 <- diff_name
rownames(diff_name2)=diff_name2[,1]
diff_name2$log10diff_namepadj <- -log10(diff_name2$FDR)
diff_name2$group <- "None"
diff_name2$group[(diff_name2$log10diff_namepadj > 1.30103) & (diff_name2$logFC > foldChange)]="Up"
diff_name2$group[(diff_name2$log10diff_namepadj > 1.30103) & (diff_name2$logFC < -foldChange)]="Down"
table(diff_name2$group)

#
a <- read.table(file = "Venn_diffname2_up&712_962交集情况.txt",sep = "\t",header = T)
my_ensembl_gene_id <- a[,1]
study_symbols<- getBM(attributes=c('ensembl_gene_id','external_gene_name',"description"), filters = 'ensembl_gene_id', values = my_ensembl_gene_id, mart = mart)
head(study_symbols)
write.csv(study_symbols,file="Venn_diffname2_up&712_962交集情况_Anno_175.csv",sep="/t")

#需要label的gene
genes = study_symbols$external_gene_name

#label genes of interest
#diff_name2$label = ""
diff_name2$label = ifelse(diff_name2$external_gene_name %in% genes, diff_name2$external_gene_name,"")

#diff_name2 to input
head(diff_name2)
input <- diff_name2[,c("ensembl_gene_id","external_gene_name","logFC","FDR","log10diff_namepadj","group","label")]
table(input$group)

#如何让你的火山图标签“乖乖听话”？个性化标签全拿捏！https://mp.weixin.qq.com/s/JpA7egexar9oN0QFjVt2EQ
#自定义颜色：
mycol <- c("#2DB2EB","#d8d8d8","#EB4232")
#自定义主题：
mytheme <- theme_classic() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        plot.margin = margin(15,5.5,5.5,5.5))

#ggplot2绘制火山图：
##foldChange=0.5849625    foldChange=log2(FC=1.5)=0.5849625
foldChange=0.5849625
padj=0.05
p <- ggplot(data = input,
            aes(x = logFC,
                y = log10diff_namepadj,
                color = group)) + #建立映射
  geom_point(size = 2.2) + #绘制散点
  scale_colour_manual(name = "", values = alpha(mycol, 0.7)) + #自定义散点颜色
  scale_x_continuous(limits = c(-15, 15),
                     breaks = seq(-15, 15, by = 3)) + #x轴限制
  scale_y_continuous(expand = expansion(add = c(2, 0)),
                     limits = c(0, 200),
                     breaks = seq(0, 200, by = 20)) + #y轴限制
  xlab("Log2FoldChange") + #X轴标签
  ylab("-Log10FDR") + #Y轴标签
  geom_hline(yintercept = c(-log10(padj)),size = 0.7,color = "black",lty = "dashed") + #水平阈值线
  geom_vline(xintercept = c(-foldChange, foldChange),size = 0.7,color = "black",lty = "dashed") + #垂直阈值线
  mytheme
p
ggsave(filename = "Part2_VolcanoPlotNewUnlabel.pdf",plot = p, width = 14,height = 12, units = "cm", dpi = 600)

label <- filter(input, label != '') %>% top_n(20, log10diff_namepadj) %>% summarise(ensembl_gene_id,
                                                                                   external_gene_name,
                                                                                   logFC,
                                                                                   FDR,
                                                                                   log10diff_namepadj,
                                                                                   group,
                                                                                   label=external_gene_name)

p5 <- p +
  geom_point(data = label,
             aes(x = logFC, y = log10diff_namepadj),
             color = '#EB4232', size = 4.5, alpha = 0.2) + #图中强调显示一下需要标注的散点
  geom_text_repel(data = label,
                  aes(x = logFC, y = log10diff_namepadj, label = label),
                  seed = 233, #可设置随机数种子
                  size = 3.5, #字体大小
                  color = 'black', #字体颜色
                  min.segment.length = 0, #始终为标签添加指引线段；若不想添加线段，则改为Inf
                  force = 2, #重叠标签间的排斥力
                  force_pull = 2, #标签和数据点间的吸引力
                  box.padding = 0.4, #标签周边填充量，默认单位为行
                  max.overlaps = Inf, #排斥重叠过多标签，设置为Inf则可以保持始终显示所有标签
                  segment.linetype = 3, #线段类型,1为实线,2-6为不同类型虚线
                  segment.color = 'black', #线段颜色
                  segment.alpha = 0.8, #线段不透明度
                  nudge_x = 15 - label$logFC, #标签x轴起始位置调整
                  direction = "y", #按y轴调整标签位置方向，若想水平对齐则为x
                  hjust = 1 #对齐标签：0右对齐，1左对齐，0.5居中
  )
p5
ggsave(filename = "Part2_VolcanoPlotNewlabel_Top20_175.pdf",plot = p5, width = 14,height = 12, units = "cm", dpi = 600)

#Heatmap 画热图
library(pheatmap)
library(RColorBrewer)
#可以用的多种颜色
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
colors <- colorRampPalette(c("navy", "white", "firebrick3"))(50)

#DEGs normalized counts table--diffExp
#diffExp=normalizeExp[diff_nameSig$ensembl_gene_id,]
head(diffExp)
input_heatmap <- diffExp[study_symbols$ensembl_gene_id,]
rownames(study_symbols) <- study_symbols[,1]
rownames(input_heatmap) <- study_symbols[rownames(input_heatmap),"external_gene_name"] 
input_heatmap_regroup = input_heatmap[,c("WT79EA", "WT79EB", "WT79EC", "NY1A", "NY2A", "NY3A", "ENY1A", "ENY2A", "ENY3A")]

write.table(input_heatmap_regroup,file="Part3_input_heatmap_regroup.xls",sep="\t",quote=F,col.names=T)         #输出差异基因校正后的表达值
nrow(input_heatmap_regroup)

heatmap=pheatmap(input_heatmap_regroup,color = colors,
                 main="175 CUT&Tag Enriched & RNA-seq Up-regulated Genes",
                 fontsize = 15,
                 scale="row",
                 border_color = NA,
                 na_col = "grey",
                 cluster_rows = T,cluster_cols = F,
                 show_rownames = T,show_colnames = T,
                 treeheight_row = 30,treeheight_col = 30,
                 cellheight = 15,cellwidth = 30,
                 cutree_row=2,cutree_col=2,
                 display_numbers = F,legend = T,
                 filename = "TPS_diffSig_Labeled.pdf"
)
heatmap


#Label基因的heatmap
##Gene_Ensemble_id
diffSig_rownames <- rownames(diffSig)
TPS_diffSig <- as.matrix(TPS[diffSig_rownames,])

##Gene_id
diffSig_rownames <- rownames(diffSig)
TPS_diffSig <- as.matrix(TPS[diffSig_rownames,])
TPS_diffSig2 <- diff_name2[match(rownames(TPS_diffSig),diff_name2$ensembl_gene_id),]
TPS_diffSig2 <- TPS_diffSig2[order(TPS_diffSig2$ensembl_gene_id),]
rownames(TPS_diffSig) <- TPS_diffSig2$external_gene_name  

##画图
heatmap=pheatmap(TPS_diffSig,color = colors,
                 main="TPS_diffSig",
                 fontsize = 15,
                 scale="row",
                 border_color = NA,
                 na_col = "grey",
                 cluster_rows = T,cluster_cols = T,
                 show_rownames = T,show_colnames = T,
                 treeheight_row = 30,treeheight_col = 30,
                 cellheight = 15,cellwidth = 30,
                 cutree_row=2,cutree_col=2,
                 display_numbers = F,legend = T,
                 filename = "TPS_diffSig.pdf"
)



#分别筛选上下调中显著性top20，作为本次测试需要添加的目标标签（共40个标签）：
up <- filter(input, group == 'Up') %>% top_n(20, log10diff_namepadj) %>% summarise(ensembl_gene_id,
                                                                                   external_gene_name,
                                                                                   logFC,
                                                                                   FDR,
                                                                                   log10diff_namepadj,
                                                                                   group,
                                                                                   label=external_gene_name)
down <- filter(input, group == 'Down') %>% top_n(20, log10diff_namepadj) %>% summarise(ensembl_gene_id,
                                                                                       external_gene_name,
                                                                                       logFC,
                                                                                       FDR,
                                                                                       log10diff_namepadj,
                                                                                       group,
                                                                                       label=external_gene_name)
head(up);head(down)

#两个函数都可用于添加标签：
##geom_text_repel():常规样式，默认参数配置
#p1 <- p +
#  geom_text_repel(data = up,
#                  aes(x = log2FoldChange, y = -log10(pvalue), label = Symbol)) +
#  geom_text_repel(data = down,
#                  aes(x = log2FoldChange, y = -log10(pvalue), label = Symbol))
#p1
##geom_label_repel():带边框样式，默认参数配置
#p2 <- p +
#  geom_label_repel(data = up,
#                   aes(x = log2FoldChange, y = -log10(pvalue), label = Symbol)) +
#  geom_label_repel(data = down,
#                   aes(x = log2FoldChange, y = -log10(pvalue), label = Symbol))
#p2

p5 <- p +
  geom_point(data = up,
             aes(x = logFC, y = log10diff_namepadj),
             color = '#EB4232', size = 4.5, alpha = 0.2) + #图中强调显示一下需要标注的散点
  geom_text_repel(data = up,
                  aes(x = logFC, y = log10diff_namepadj, label = label),
                  seed = 233, #可设置随机数种子
                  size = 3.5, #字体大小
                  color = 'black', #字体颜色
                  min.segment.length = 0, #始终为标签添加指引线段；若不想添加线段，则改为Inf
                  force = 2, #重叠标签间的排斥力
                  force_pull = 2, #标签和数据点间的吸引力
                  box.padding = 0.4, #标签周边填充量，默认单位为行
                  max.overlaps = Inf, #排斥重叠过多标签，设置为Inf则可以保持始终显示所有标签
                  segment.linetype = 3, #线段类型,1为实线,2-6为不同类型虚线
                  segment.color = 'black', #线段颜色
                  segment.alpha = 0.8, #线段不透明度
                  nudge_x = 9 - up$logFC, #标签x轴起始位置调整
                  direction = "y", #按y轴调整标签位置方向，若想水平对齐则为x
                  hjust = 1 #对齐标签：0右对齐，1左对齐，0.5居中
  )
p5

#左侧添加方法相同：
p6 <- p5 +
  geom_point(data = down,
             aes(x = logFC, y = log10diff_namepadj),color = '#2DB2EB',size = 4.5, alpha = 0.2) +
  geom_text_repel(data = down,
                  aes(x = logFC, y = log10diff_namepadj, label = label),
                  seed = 233,
                  size = 3.5,
                  color = 'black',
                  min.segment.length = 0,
                  force = 2,
                  force_pull = 2,
                  box.padding = 0.4,
                  max.overlaps = Inf,
                  segment.linetype = 3,
                  segment.color = 'black',
                  segment.alpha = 0.8,
                  nudge_x = -8 - down$logFC,
                  direction = "y",
                  hjust = 1 #改为左对齐即可
  )
p6
ggsave(filename = "Part2_VolcanoPlotNew.pdf",plot = p6, width = 14,height = 12, units = "cm", dpi = 600)


# Step 3 差异表达基因功能富集分析
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(clusterProfiler)
library(DOSE)
library(ggplot2)
library(stringr)
library(AnnotationDbi)
library(org.Dm.eg.db)
library(Cairo)
library(enrichplot)
#GO Enrichment 
head(diff_nameSig)
#keyType = "ENSEMBL", 这里只支持ENSEMBL
ALL <- enrichGO(gene = diff_nameSig$ensembl_gene_id, 
                OrgDb = org.Dm.eg.db, 
                keyType = "ENSEMBL",
                ont = 'ALL',
                pvalueCutoff  = 0.05,
                pAdjustMethod = "BH",  
                qvalueCutoff  = 0.1, readable=T)  #一步到位
				
write.csv(as.data.frame(ALL@result), file="Part2_GOALL.csv",sep="\t", quote=FALSE)
GOFile = read.csv("Part2_GOALL.csv", header = T, sep="\t")

Goplot <- dotplot(ALL, split = "ONTOLOGY", font.size = 8, showCategory = 5) + 
				facet_grid(ONTOLOGY ~ ., scale = "free") + 
				scale_y_discrete(labels = function(x) str_wrap(x, width =50)) + 
                scale_size(range=c(2, 6))   #设置点的大小

ggsave("Part2_GOplot.pdf", plot = Goplot, width = 8, height = 8) 
head(ALL,1);dim(ALL)				
sum(ALL$ONTOLOGY=="BP") #Biological process基因产物参与的生物路径或机制
sum(ALL$ONTOLOGY=="CC") #Cellular component基因产物在细胞内外的位置
sum(ALL$ONTOLOGY=="MF") #Molecular function基因产物分子层次的功能

################################################################################error in kegg
diff_nameSig <- read.table(file = "Part1_diff_nameSig.xls",sep = "\t",header=T)
library(stringr)
library(DOSE)
data(geneList, package="DOSE")
columns(org.Dm.eg.db)
gene.df <- bitr(diff_nameSig$ensembl_gene_id, fromType = "ENSEMBL", 
                toType = c("SYMBOL","ENTREZID"),
                OrgDb = org.Dm.eg.db) 

gene.kegg <- bitr_kegg(gene.df$ENTREZID,fromType="ncbi-geneid",
                       toType="kegg",organism='dme')
head(gene.kegg)


getOption("clusterProfiler.download.method")

utils::download.file(url, quiet = quiet, method = "libcurl")
options(clusterProfiler.download.method = "auto")
#options(clusterProfiler.download.method = "libcurl")
#options(clusterProfiler.download.method = "wininet")

ekegg <- enrichKEGG(gene = gene.kegg$kegg, 
                    organism = "dme", 
                    keyType = "kegg",
                    pvalueCutoff  = 0.05,
                    pAdjustMethod = "BH",  
                    qvalueCutoff  = 0.1)  #一步到位

erich.go.BP = enrichGO(gene = DEG.gene_symbol,
                       OrgDb = org.Dm.eg.db,
                       keyType = "SYMBOL",
                       ont = "BP",
                       pvalueCutoff = 0.5,
                       qvalueCutoff = 0.5)

erich.go.CC = enrichGO(gene = DEG.gene_symbol,
                       OrgDb = org.Dm.eg.db,
                       keyType = "SYMBOL",
                       ont = "CC",
                       pvalueCutoff = 0.5,
                       qvalueCutoff = 0.5)

erich.go.MF = enrichGO(gene = DEG.gene_symbol,
                       OrgDb = org.Dm.eg.db,
                       keyType = "SYMBOL",
                       ont = "MF",
                       pvalueCutoff = 0.5,
                       qvalueCutoff = 0.5)

dotplot(erich.go.CC)
dotplot(erich.go.BP)
dotplot(erich.go.MF)

# save image to file
pdf(file="./20200926-enrich.go.BP.Dotplot.pdf",width = 10,height = 6	)
dotplot(erich.go.BP)
dev.off()

# KEGG analysis
#查看数据库中的所有名字
columns(org.Dm.eg.db)

# convert id
DEG.gene_symbol$DEG.entrez_id = mapIds(x = org.Dm.eg.db,
                       keys = DEG.gene_symbol$DEG.gene_symbol,
                       keytype = "SYMBOL",
                       column = "ENTREZID")

DEG.gene_symbol$DEG.uniprot_id = mapIds(x = org.Dm.eg.db,
                       keys = DEG.gene_symbol$DEG.gene_symbol,
                       keytype = "SYMBOL",
                       column = "UNIPROT")

erich.kegg.res <- enrichKEGG(gene = DEG.gene_symbol$DEG.entrez_id,
                             organism = "dme",
                             keyType = "kegg",
                             pAdjustMethod = "BH")

barplot(erich.kegg.res)








