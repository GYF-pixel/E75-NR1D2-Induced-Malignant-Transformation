#1.FastQC
cd /storage/maxianjueLab/guoyifan/CUT_Tag

projPath="/storage/maxianjueLab/guoyifan/CUT_Tag"
module load fastqc/0.11.9

#!/bin/bash
#SBATCH -J fastqc
#SBATCH -p amd-ep2
#SBATCH -q normal
#SBATCH --mem-per-cpu=358400
sh /storage/maxianjueLab/guoyifan/CUT_Tag/fastqFileQC/fastqc.sh 
fastqc -o ${projPath}/fastqFileQC ${projPath}/Rawdata/E75/*.fq.gz

#!/bin/bash
#SBATCH -J bowtie2
#SBATCH -p amd-ep2
#SBATCH -q normal
#SBATCH --mem=358400

cd /storage/maxianjueLab/guoyifan/CUT_Tag

cd /storage/maxianjueLab/guoyifan/CUT_Tag/Rawdata/E75

module load bowtie/2.4.2
projPath="/storage/maxianjueLab/guoyifan/CUT_Tag"
 mkdir  -p ${projPath}/alignment/sam/bowtie2_summary
 mkdir  -p ${projPath}/alignment/bam
 mkdir  -p ${projPath}/alignment/bed
 mkdir  -p ${projPath}/alignment/bedgraph

#2. Bowtie2 Alignment to Drosophila_melanogaster & Ecoli &DNA
##== linux 命令 ==##

module load samtools/1.11
module load bowtie/2.4.2
projPath="/storage/maxianjueLab/guoyifan/CUT_Tag"
#Drosophila_melanogaster
ref="/storage/maxianjueLab/guoyifan/species_reference/DroMelanogaster/Drosophila_melanogaster.BDGP6.32"
#Ecoli
spikeInRef="/storage/maxianjueLab/guoyifan/species_reference/Escherichia_coli/Ecoli"
#DNA standard sample
spikeInRef="/storage/maxianjueLab/guoyifan/species_reference/CUTTagSpikeInDNA/SpikeInDNA"
cores=24
## Build the bowtie2 reference genome index if needed:
## bowtie2-build path/to/hg38/fasta/hg38.fa /path/to/bowtie2Index/hg38
## bowtie2-build /storage/maxianjueLab/guoyifan/species_reference/Drosophila_melanogaster_BDGP6_32/Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa /storage/maxianjueLab/guoyifan/species_reference/DroMelanogaster/Drosophila_melanogaster.BDGP6.32
## bowtie2-build /storage/maxianjueLab/guoyifan/species_reference/Escherichia_coli/GCF_000005845.2_ASM584v2_genomic.fna /storage/maxianjueLab/guoyifan/species_reference/Escherichia_coli/Ecoli
## bowtie2-build /storage/maxianjueLab/guoyifan/species_reference/CUTTagSpikeInDNA/CUTTagSpikeInDNA.fna /storage/maxianjueLab/guoyifan/species_reference/CUTTagSpikeInDNA/SpikeInDNA


for histName in E75_Flag1 E75_Flag2 E75_Flag3 E75_Flag4 ENY_N10 ENY_N2 ENY_N6 ENY_N8 ENY_Y2 ENY_Y4 ENY_Y6 ENY_Y8 NY_N1 NY_N5 NY_N7 NY_N9 NY_Y1 NY_Y3 NY_Y5 NY_Y7;do
  {
bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 200 -X 800 -p ${cores} -x ${ref} \
-1 ${projPath}/00Rawdata/${histName}_R1.fq.gz \
-2 ${projPath}/00Rawdata/${histName}_R2.fq.gz \
-S ${projPath}/alignment/sam/${histName}_bowtie2.sam &> ${projPath}/alignment/sam/bowtie2_summary/${histName}_bowtie2.txt

bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 200 -X 800 -p ${cores} -x ${spikeInRef} \
-1 ${projPath}/00Rawdata/${histName}_R1.fq.gz \
-2 ${projPath}/00Rawdata/${histName}_R2.fq.gz \
-S ${projPath}/alignment/sam/${histName}_bowtie2_spikeIn.sam &> ${projPath}/alignment/sam/bowtie2_summary/${histName}_bowtie2_spikeIn.txt

  seqDepthDouble=`samtools view -F 0x4 $projPath/alignment/sam/${histName}_bowtie2_spikeIn.sam | wc -l`
  seqDepth=$((seqDepthDouble/2))
  echo $seqDepth >$projPath/alignment/sam/bowtie2_summary/${histName}_bowtie2_spikeIn.seqDepth
  }&
done

## -I 是最短序列  -X 是最长序列 -p 线程 -x 指定参考基因组的index -1 序列1文件 -2 序列2文件 -S 生成的sam文件
##由于文件比较多，可以写一个简单的for循环的shell脚本放到文件目录下面（ for i in `cat file_list` ; do done）
##比对到研究物种的参考基因组上，其中*-rep*是加了研究抗体的，IgG_rep1是IgG，用来call峰的时候去除背景的，因为目前只做了一个样本，4个重复，因此后面是没有办法做差异peak的，所以在后面用基因注释的时候选用的是全部的peak做的基因注释。


#3. SAM to BAM to BED
module load samtools/1.11
module load bedtools/2.30.0
module load bowtie/2.4.2
projPath="/storage/maxianjueLab/guoyifan/CUT_Tag"
## Filter and keep the mapped read pairs
## 筛选和保留比对上的双端 reads 
for i in E75_Flag1 E75_Flag2 E75_Flag3 E75_Flag4 ENY_N10 ENY_N2 ENY_N6 ENY_N8 ENY_Y2 ENY_Y4 ENY_Y6 ENY_Y8 NY_N1 NY_N5 NY_N7 NY_N9 NY_Y1 NY_Y3 NY_Y5 NY_Y7;do
  {
  samtools view -bS -F 0x4 $projPath/alignment/sam/${i}_bowtie2.sam > $projPath/alignment/bam/${i}_bowtie2.mapped.bam
## Convert into bed file format
## 将 BAM 文件转换为 bed 文件格式
  bedtools bamtobed -i $projPath/alignment/bam/${i}_bowtie2.mapped.bam -bedpe > $projPath/alignment/bed/${i}_bowtie2.bed
## Keep the read pairs that are on the same chromosome and fragment length less than 1000bp.
## 保留那些在同一条染色体且片段长度小于 1000bp 的双端 reads
  awk '$1==$4 && $6-$2 < 1000 {print $0}' $projPath/alignment/bed/${i}_bowtie2.bed > $projPath/alignment/bed/${i}_bowtie2.clean.bed
## Only extract the fragment related columns
## 仅提取片段相关的列 
  cut -f 1,2,6 $projPath/alignment/bed/${i}_bowtie2.clean.bed | sort -k1,1 -k2,2n -k3,3n  > $projPath/alignment/bed/${i}_bowtie2.fragments.bed
  }&
done

##3.1 评估重复性
##为了研究重复之间和不同条件下的可重复性，基因组被分成 500 bp bin，每个 bin reads 计数的 log2 转换值的皮尔逊相关性在重复数据集之间计算。多个重复和 IgG 对照数据集显示在层次化聚类关联矩阵中。
for i in E75_Flag1 E75_Flag2 E75_Flag3 E75_Flag4 ENY_N10 ENY_N2 ENY_N6 ENY_N8 ENY_Y2 ENY_Y4 ENY_Y6 ENY_Y8 NY_N1 NY_N5 NY_N7 NY_N9 NY_Y1 NY_Y3 NY_Y5 NY_Y7; do
  {
  binLen=500
  awk -v w=$binLen '{print $1, int(($2 + $3)/(2*w))*w + w/2}' $projPath/alignment/bed/${i}_bowtie2.fragments.bed |\
  sort -k1,1V -k2,2n |\
  uniq -c |\
  awk -v OFS="\t" '{print $2, $3, $1}' |\
  sort -k1,1V -k2,2n  > $projPath/alignment/bed/${i}_bowtie2.fragmentsCount.bin$binLen.bed
  }&
done


##3.2 标准化
module load samtools/1.11
module load bedtools/2.30.0
module load bowtie/2.4.2
module load R/4.2.1
projPath="/storage/maxianjueLab/guoyifan/CUT_Tag"
chromSize="/storage/maxianjueLab/guoyifan/species_reference/Drosophila_melanogaster_BDGP6_32/Drosophila_melanogaster.chrom.size"
#3.2.1 生成chromsize文件
#samtools faidx ref.genome ##后面的bedtools -g的内容通过samtools的faidx来进行构建获得索引文件
#cd /storage/maxianjueLab/guoyifan/species_reference/Drosophila_melanogaster_BDGP6_32
#samtools faidx Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa
#cut -f1,2 Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa.fai > Drosophila_melanogaster.chrom.size
#chromSize="/storage/maxianjueLab/guoyifan/species_reference/Drosophila_melanogaster_BDGP6_32/Drosophila_melanogaster.chrom.size"

#3.2.2 标准化得到bedgraph文件放到IGV里查看
##这里的seqDepth是前面的数据比对处理当中获得测序深度值，可以调取txt文件，来进行后续的循环处理
## -lt less than -gt greater than
## `` backtick  A backtick is not a quotation sign. It has a very special meaning. Everything you type between backticks is evaluated (executed) by the shell before the main command https://unix.stackexchange.com/questions/27428/what-does-backquote-backtick-mean-in-commands
## " " quotation sign

for histName in E75_Flag1 E75_Flag2 E75_Flag3 E75_Flag4 ENY_N10 ENY_N2 ENY_N6 ENY_N8 ENY_Y2 ENY_Y4 ENY_Y6 ENY_Y8 NY_N1 NY_N5 NY_N7 NY_N9 NY_Y1 NY_Y3 NY_Y5 NY_Y7; do
  {
seqDepthDouble=`samtools view -F 0x4 $projPath/alignment/sam/${histName}_bowtie2_spikeIn.sam | wc -l`
seqDepth=$((seqDepthDouble/2))
echo $seqDepth >$projPath/alignment/sam/bowtie2_summary/${histName}_bowtie2_spikeIn.seqDepth
  if [[ "$seqDepth" -gt "1" ]]; then
    scale_factor=`echo "10000 / $seqDepth" | bc -l`
    echo "Scaling factor for $histName is: $scale_factor!"
    bedtools genomecov -bg -scale $scale_factor -i $projPath/alignment/bed/${histName}_bowtie2.fragments.bed -g $chromSize > $projPath/alignment/bedgraph/${histName}_bowtie2.fragments.normalized.bedgraph  
  fi  
  }&
done


#4. 画热图Heatmap visualization on specific regions
##== linux command ==##

#mkdir -p $projPath/alignment/bigwig     
module load samtools/1.11
module load bedtools/2.30.0
module load bowtie/2.4.2
module load R/4.2.1
source ~/.bashrc
conda activate deeptools
projPath="/storage/maxianjueLab/guoyifan/CUT_Tag"
cd /storage/maxianjueLab/guoyifan/CUT_Tag/alignment/bam

#4.1 先对bam文件转换成bw文件 https://www.jianshu.com/p/491557586118
 for i in E75_Flag1 E75_Flag2 E75_Flag3 E75_Flag4 ENY_N10 ENY_N2 ENY_N6 ENY_N8 ENY_Y2 ENY_Y4 ENY_Y6 ENY_Y8 NY_N1 NY_N5 NY_N7 NY_N9 NY_Y1 NY_Y3 NY_Y5 NY_Y7;do
   {
     samtools sort -@ 12 $projPath/alignment/bam/${i}_bowtie2.mapped.bam -o $projPath/alignment/bam/${i}_bowtie2.mapped.sorted.bam
     samtools index $projPath/alignment/bam/${i}_bowtie2.mapped.sorted.bam
	 bamCoverage -b $projPath/alignment/bam/${i}_bowtie2.mapped.sorted.bam -o $projPath/alignment/bigwig/${i}_raw.bw
	 bamCoverage -b $projPath/alignment/bam/${i}_bowtie2.mapped.sorted.bam -o $projPath/alignment/bigwig/${i}_normalized_BPM_50bin.bw --binSize 50 --normalizeUsing BPM	 
}&
done

#3.2 Picard remove duplication
#PCR cycle == 13, Picard procedure is omitted
#java -Xms40g -Xmx40g -XX:ParallelGCThreads=12 -jar {PICARD} MarkDuplicates I={input} O={output[0]} M={output[1]} ASO=coordinate REMOVE_DUPLICATES=true 2>{log}

##3.3 Sort & MAPQ20
 for i in E75_Flag1 E75_Flag2 E75_Flag3 E75_Flag4 ENY_N10 ENY_N2 ENY_N6 ENY_N8 ENY_Y2 ENY_Y4 ENY_Y6 ENY_Y8 NY_N1 NY_N5 NY_N7 NY_N9 NY_Y1 NY_Y3 NY_Y5 NY_Y7;do
   {
     samtools sort -@ 12 $projPath/alignment/bam/${i}_bowtie2.mapped.bam -o $projPath/alignment/bam/${i}_bowtie2.mapped.sorted.bam
	 samtools view -t 12 -h -f bam -F "mapping_quality >= 20" $projPath/alignment/bam/${i}_bowtie2.mapped.bam -o $projPath/alignment/bam/${i}_bowtie2.mapped.sorted.mapq.bam
	 samtools index $projPath/alignment/bam/${i}_bowtie2.mapped.sorted.mapq.bam
	 }&
done

#3.4 bam to bw
 for i in E75_Flag1 E75_Flag2 E75_Flag3 E75_Flag4 ENY_N10 ENY_N2 ENY_N6 ENY_N8 ENY_Y2 ENY_Y4 ENY_Y6 ENY_Y8 NY_N1 NY_N5 NY_N7 NY_N9 NY_Y1 NY_Y3 NY_Y5 NY_Y7;do
   {
	 bamCoverage -b $projPath/alignment/bam/${i}_bowtie2.mapped.sorted.mapq.bam -o $projPath/alignment/bigwig/${i}_raw.bw
	 bamCoverage -b $projPath/alignment/bam/${i}_bowtie2.mapped.sorted.mapq.bam -o $projPath/alignment/bigwig/${i}_normalized_BPM_50bin.bw --binSize 50 --normalizeUsing BPM	 
}&
done

#RPKM = Reads Per Kilobase per Million mapped reads;
#CPM= Counts Per Million mapped reads, same as CPM in RNA-seq;
#BPM = Bins Per Million mapped reads, same as TPM in RNA-seq;
#RPGC = reads per genomic content (1x normalization); Mapped reads are considered after blacklist filtering (if applied).
#需要注意的是选择了RPGC后，需要提供--effectiveGenomeSize

#PCR循环是13，所以没用Picard来去重

#MAPQ20
 for i in E75_Flag1 E75_Flag2 E75_Flag3 E75_Flag4 ENY_N10 ENY_N2 ENY_N6 ENY_N8 ENY_Y2 ENY_Y4 ENY_Y6 ENY_Y8 NY_N1 NY_N5 NY_N7 NY_N9 NY_Y1 NY_Y3 NY_Y5 NY_Y7;do
   {
	 samtools view -t 12 -h -f bam -F "mapping_quality >= 20" -o {output} {input}
   	 }&
done

rule MAPQ20_select_bam:
    input:
        "bam.bt2/{sample}_bt2_hg38_sort_rmdup.bam"
    output:
        "bam.bt2/{sample}_bt2_hg38_sort_rmdup_MAPQ20.bam"
    shell:
        """sambamba view -t 12 -h -f bam -F "mapping_quality >= 20" -o {output} {input}"""

## 筛选和保留比对上的双端 reads 
for i in E75_Flag1 E75_Flag2 E75_Flag3 E75_Flag4 ENY_N10 ENY_N2 ENY_N6 ENY_N8 ENY_Y2 ENY_Y4 ENY_Y6 ENY_Y8 NY_N1 NY_N5 NY_N7 NY_N9 NY_Y1 NY_Y3 NY_Y5 NY_Y7;do
  {
  samtools view -bS -F 0x4 $projPath/alignment/sam/${i}_bowtie2.sam > $projPath/alignment/bam/${i}_bowtie2.mapped.bam
## Convert into bed file format
## 将 BAM 文件转换为 bed 文件格式
  bedtools bamtobed -i $projPath/alignment/bam/${i}_bowtie2.mapped.bam -bedpe > $projPath/alignment/bed/${i}_bowtie2.bed
## Keep the read pairs that are on the same chromosome and fragment length less than 1000bp.
## 保留那些在同一条染色体且片段长度小于 1000bp 的双端 reads
  awk '$1==$4 && $6-$2 < 1000 {print $0}' $projPath/alignment/bed/${i}_bowtie2.bed > $projPath/alignment/bed/${i}_bowtie2.clean.bed
## Only extract the fragment related columns
## 仅提取片段相关的列 
  cut -f 1,2,6 $projPath/alignment/bed/${i}_bowtie2.clean.bed | sort -k1,1 -k2,2n -k3,3n  > $projPath/alignment/bed/${i}_bowtie2.fragments.bed
  }&
done

# align stats 
rule align_stats_final:
    input:
        "bam.bt2/{sample}_bt2_hg38_sort_rmdup_MAPQ20.bam"
    output:
        "bam.bt2/{sample}_bt2_hg38_sort_rmdup_MAPQ20.mapping_stats"
    shell:
        "samtools stats -@ 5 {input} > {output}"




#4.2 对注释文件进行提取来获得基因的注释内容
###文件重命名
#for histName in E75_Flag1 E75_Flag2 E75_Flag4 ENY_N1 ENY_N2 ENY_N3 ENY_Y1 ENY_Y2 ENY_Y5 NY_N1 NY_N2 NY_N3 NY_Y1 NY_Y3 NY_Y4; do
#{
#rename ${histName}_refPoint_TSS_heatmap ${histName}_refPoint_TSS_TES_heatmap ${histName}_refPoint_TSS_heatmap.pdf
#rename ${histName}_refPoint_TSS_data ${histName}_refPoint_TSS_TES_data ${histName}_refPoint_TSS_data.gz
#   }&
#done


##reference-point TSS
##AD002A 暗红色
#cd /storage/maxianjueLab/guoyifan/species_reference/Drosophila_melanogaster_BDGP6_32
#cat Drosophila_melanogaster.BDGP6.32.109.gtf | cut -f1,4,5,9 | cut -f1 -d";" | awk '{print $1, $2, $3, $5}' | sed -e 's/ /\t/g' | sed -e 's/\"//g' | sed -e 's/transcript\://g' > Drosophila_melanogaster_BDGP6_32_deeptools.genomic.bed

ref_gtf="/storage/maxianjueLab/guoyifan/species_reference/Drosophila_melanogaster_BDGP6_32/Drosophila_melanogaster.BDGP6.32.109.gtf"
projPath="/storage/maxianjueLab/guoyifan/CUT_Tag"

for histName in ENY_N10 ENY_N2 ENY_N6 ENY_N8 ENY_Y2 ENY_Y4 ENY_Y6 ENY_Y8 NY_N1 NY_N5 NY_N7 NY_N9 NY_Y1 NY_Y3 NY_Y5 NY_Y7; do
{
computeMatrix reference-point -S $projPath/alignment/bigwig/${histName}_normalized_BPM_50bin.bw \
							   -p 24 \
							  --binSize 50 \
							  --regionBodyLength 5000 \
							  --referencePoint TSS \
							  --afterRegionStartLength 3000 \
							  --beforeRegionStartLength 3000 \
							  -R $ref_gtf \
							  --skipZeros \
							  -o $projPath/deeptools/${histName}_refPoint_TSS_TES_data.gz

plotHeatmap -m $projPath/deeptools/${histName}_refPoint_TSS_TES_data.gz \
            --missingDataColor 1 \
            --colorList 'white,#AD002A' \
            --heatmapHeight 12 \
			--sortUsing sum --startLabel "TSS" \
            --endLabel "TES" --xAxisLabel "" \
            --regionsLabel "Peaks" \
            --samplesLabel "${histName}" \
            -o $projPath/deeptools/${histName}_refPoint_TSS_TES_heatmap.pdf
    }&
done



for histName in E75_Flag1 E75_Flag2 E75_Flag3 E75_Flag4; do
{
computeMatrix reference-point -S $projPath/alignment/bigwig/${histName}_normalized_BPM_50bin.bw \
							   -p 24 \
							  --binSize 50 \
							  --regionBodyLength 5000 \
							  --referencePoint TSS \
							  --afterRegionStartLength 3000 \
							  --beforeRegionStartLength 3000 \
							  -R $ref_gtf \
							  --skipZeros \
							  -o $projPath/deeptools/${histName}_refPoint_TSS_TES_data.gz

plotHeatmap -m $projPath/deeptools/${histName}_refPoint_TSS_TES_data.gz \
            --missingDataColor 1 \
            --colorList 'white,#AD002A' \
            --heatmapHeight 12 \
			--sortUsing sum --startLabel "TSS" \
            --endLabel "TES" --xAxisLabel "" \
            --regionsLabel "Peaks" \
            --samplesLabel "${histName}" \
            -o $projPath/deeptools/${histName}_refPoint_TSS_TES_heatmap.pdf
    }&
done




###925E9F 暗紫色
ref_gtf="/storage/maxianjueLab/guoyifan/species_reference/Drosophila_melanogaster_BDGP6_32/Drosophila_melanogaster.BDGP6.32.109.gtf"
projPath="/storage/maxianjueLab/guoyifan/CUT_Tag"
for histName in E75_Flag1 E75_Flag2 E75_Flag4 ENY_N1 ENY_N2 ENY_N3 ENY_Y1 ENY_Y2 ENY_Y5 NY_N1 NY_N2 NY_N3 NY_Y1 NY_Y3 NY_Y4; do
{
computeMatrix reference-point -S $projPath/alignment/bigwig/${histName}_normalized_BPM_50bin.bw \
							   -p 24 \
							  --referencePoint TSS \
							  --afterRegionStartLength 3000 \
							  --beforeRegionStartLength 3000 \
							  -R $ref_gtf \
							  --skipZeros \
							  -o $projPath/deeptools/${histName}_refPoint_TSS_data.gz

plotHeatmap -m $projPath/deeptools/${histName}_refPoint_TSS_data.gz \
            --missingDataColor 1 \
            --colorList 'white,#925E9F' \
            --heatmapHeight 12 \
			--sortUsing sum --startLabel "TSS" \
            --endLabel "TES" --xAxisLabel "" \
            --regionsLabel "Peaks" \
            --samplesLabel "${histName}" \
            -o $projPath/deeptools/${histName}_refPoint_TSS_heatmap.pdf
    }&
done

for histName in E75 ENY_N ENY_Y NY_N NY_Y; do
{
computeMatrix reference-point -S $projPath/alignment/bigwig/${histName}_normalized_BPM_50bin.bw \
							   -p 6 \
							  --referencePoint TSS \
							  --afterRegionStartLength 3000 \
							  --beforeRegionStartLength 3000 \
							  -R $ref_gtf \
							  --skipZeros \
							  -o $projPath/deeptools/${histName}_refPoint_TSS_data.gz

plotHeatmap -m $projPath/deeptools/${histName}_refPoint_TSS_data.gz \
            --missingDataColor 1 \
            --colorList 'white,#925E9F' \
            --heatmapHeight 12 \
			--sortUsing sum --startLabel "TSS" \
            --endLabel "TES" --xAxisLabel "" \
            --regionsLabel "Peaks" \
            --samplesLabel "${histName}" \
            -o $projPath/deeptools/${histName}_refPoint_TSS_heatmap.pdf
    }&
done


##scale-regions
ref_from_gtf="/storage/maxianjueLab/guoyifan/species_reference/Drosophila_melanogaster_BDGP6_32/Drosophila_melanogaster.BDGP6.32.109.gtf"
projPath="/storage/maxianjueLab/guoyifan/CUT_Tag"
histName="E75"
computeMatrix scale-regions -S $projPath/alignment/bigwig/E75_Flag1_raw.bw \
                               $projPath/alignment/bigwig/E75_Flag2_raw.bw \
                               $projPath/alignment/bigwig/E75_Flag4_raw.bw \
							   -p 10 \
							  --binSize 50 \
							  --regionBodyLength 5000 \
							  --afterRegionStartLength 3000 \
							  --beforeRegionStartLength 3000 \
							  -R $ref_from_gtf \
							  --skipZeros \
							  -o $projPath/deeptools/${histName}_combine-scaleRegion-data.gz
							  
plotHeatmap -m $projPath/deeptools/${histName}_combine-scaleRegion-data.gz \
            --missingDataColor 1 \
            --colorList 'white,#0066CC' \
            --heatmapHeight 12 \
            -o $projPath/deeptools/${histName}_scaleRegion_heatmap.pdf


#5.1 Peak Calling MACS2
##== linux 命令 ==##
#mkdir -p $projPath/peakCalling/MACS2 https://mp.weixin.qq.com/s/IjbRisuSnfUWsYs1dSn5_g
projPath="/storage/maxianjueLab/guoyifan/CUT_Tag"
module load samtools/1.11
module load bedtools/2.30.0
module load bowtie/2.4.2
module load R/4.2.1
source ~/.bashrc
conda activate macs2
# No control
for histName in E75_Flag1 E75_Flag2 E75_Flag4 ENY_N1 ENY_N2 ENY_N3 ENY_Y1 ENY_Y2 ENY_Y5 NY_N1 NY_N2 NY_N3 NY_Y1 NY_Y3 NY_Y4; do
  {
macs2 callpeak -t ${projPath}/alignment/bam/${histName}_bowtie2.mapped.sorted.bam \
      -g dm -f BAMPE -n macs2_${histName}_peak_q0.1 --outdir $projPath/peakCalling/MACS2 -q 0.1 --keep-dup all 2>${projPath}/peakCalling/MACS2/macs2Peak_${histName}_summary.txt
    }&
done

for histName in E75_Flag1 E75_Flag2 E75_Flag3 E75_Flag4 ENY_N10 ENY_N2 ENY_N6 ENY_N8 ENY_Y2 ENY_Y4 ENY_Y6 ENY_Y8 NY_N1 NY_N5 NY_N7 NY_N9 NY_Y1 NY_Y3 NY_Y5 NY_Y7; do
  {
macs2 callpeak -t ${projPath}/alignment/bam/${histName}_bowtie2.mapped.sorted.bam \
      -g dm -f BAMPE -n macs2_${histName}_peak_q0.1 --outdir $projPath/peakCalling/MACS2 -q 0.1 --keep-dup all 2>${projPath}/peakCalling/MACS2/macs2Peak_${histName}_summary.txt
    }&
done

# With control
for histName in E75_Flag1 E75_Flag2 E75_Flag4 ENY_N1 ENY_N2 ENY_N3 ENY_Y1 ENY_Y2 ENY_Y5 NY_N1 NY_N2 NY_N3 NY_Y1 NY_Y3 NY_Y4; do
  {
	controlName="CT"
		macs2 callpeak -t ${projPath}/alignment/bam/${histName}_bowtie2.mapped.sorted.bam \
		-c ${projPath}/alignment/bam/${controlName}_bowtie2.mapped.sorted.bam \
		-g dm -f BAMPE -n macs2_${histName}_peak_q0.1 --outdir $projPath/peakCalling/MACS2 -q 0.1 --keep-dup all 2>${projPath}/peakCalling/MACS2/macs2Peak_${histName}_summary.txt
    }&
done



 #5.2 Peak calling SEACR

##== linux 命令==##
##Peak calling
#mkdir -p $projPath/peakCalling/SEACR
module load samtools/1.11
module load bedtools/2.30.0
module load bowtie/2.4.2
module load R/4.2.1
##== linux 命令==##
projPath="/storage/maxianjueLab/guoyifan/CUT_Tag"
seacr="/storage/maxianjueLab/guoyifan/Software/SEACR-master/SEACR_1.3.sh"
histControl="CT1"
#With Control
for histName in E75_Flag1 E75_Flag2 E75_Flag4 ENY_N1 ENY_N2 ENY_N3 ENY_Y1 ENY_Y2 ENY_Y5 NY_N1 NY_N2 NY_N3 NY_Y1 NY_Y3 NY_Y4; do
  {
  bash $seacr $projPath/alignment/bedgraph/${histName}_bowtie2.fragments.normalized.bedgraph \
            $projPath/alignment/bedgraph/${histControl}_bowtie2.fragments.normalized.bedgraph \
            non stringent $projPath/peakCalling/SEACR/${histName}_seacr_control.peaks
  bash $seacr $projPath/alignment/bedgraph/${histName}_bowtie2.fragments.normalized.bedgraph 0.01 non stringent $projPath/peakCalling/SEACR/${histName}_seacr_top0.01.peaks
  bash $seacr $projPath/alignment/bedgraph/${histName}_bowtie2.fragments.normalized.bedgraph 0.05 non stringent $projPath/peakCalling/SEACR/${histName}_seacr_top0.05.peaks
    }&
done


#Without Control
for histName in E75_Flag1 E75_Flag2 E75_Flag4 ENY_N1 ENY_N2 ENY_N3 ENY_Y1 ENY_Y2 ENY_Y5 NY_N1 NY_N2 NY_N3 NY_Y1 NY_Y3 NY_Y4; do
  {
  bash $seacr $projPath/alignment/bedgraph/${histName}_bowtie2.fragments.normalized.bedgraph 0.01 non stringent $projPath/peakCalling/SEACR/${histName}_seacr_top0.01.peaks
  bash $seacr $projPath/alignment/bedgraph/${histName}_bowtie2.fragments.normalized.bedgraph 0.05 non stringent $projPath/peakCalling/SEACR/${histName}_seacr_top0.05.peaks
    }&
done



##referencePoint center SEACR  
#339933 暗绿色
ref_gtf="/storage/maxianjueLab/guoyifan/species_reference/Drosophila_melanogaster_BDGP6_32/Drosophila_melanogaster.BDGP6.32.109.gtf"
projPath="/storage/maxianjueLab/guoyifan/CUT_Tag"
for histName in E75_Flag1 E75_Flag2 E75_Flag4 ENY_N1 ENY_N2 ENY_N3 ENY_Y1 ENY_Y2 ENY_Y5 NY_N1 NY_N2 NY_N3 NY_Y1 NY_Y3 NY_Y4; do
{
awk '{split($6, summit, ":"); split(summit[2], region, "-"); print summit[1]"\t"region[1]"\t"region[2]}' $projPath/peakCalling/SEACR/${histName}_seacr_top0.01.peaks.stringent.bed >$projPath/peakCalling/SEACR/${histName}_seacr_top0.01.peaks.summitRegion.bed
computeMatrix reference-point -S $projPath/alignment/bigwig/${histName}_normalized_BPM_50bin.bw \
							   -p 24 \
							  --referencePoint center \
							  --afterRegionStartLength 3000 \
							  --beforeRegionStartLength 3000 \
							  -R $projPath/peakCalling/SEACR/${histName}_seacr_top0.01.peaks.summitRegion.bed \
							  --skipZeros \
							  -o $projPath/deeptools/${histName}_SEACR_refPoint_center_data.gz

plotHeatmap -m $projPath/deeptools/${histName}_SEACR_refPoint_center_data.gz \
            --missingDataColor 1 \
            --colorList 'white,#339933' \
            --heatmapHeight 12 \
			--sortUsing sum --startLabel "Peak Start" \
            --endLabel "Peak End" --xAxisLabel "" \
            --regionsLabel "Peaks" \
            --samplesLabel "${histName}" \
            -o $projPath/deeptools/${histName}_SEACR_refPoint_center_heatmap.pdf
  }&
done



##referencePoint center MACS2 summits
#0066CC 浅蓝色
ref_gtf="/storage/maxianjueLab/guoyifan/species_reference/Drosophila_melanogaster_BDGP6_32/Drosophila_melanogaster.BDGP6.32.109.gtf"
projPath="/storage/maxianjueLab/guoyifan/CUT_Tag"
for histName in E75_Flag1 E75_Flag2 E75_Flag4 ENY_N1 ENY_N2 ENY_N3 ENY_Y1 ENY_Y2 ENY_Y5 NY_N1 NY_N2 NY_N3 NY_Y1 NY_Y3 NY_Y4; do
{
#awk '{split($6, summit, ":"); split(summit[2], region, "-"); print summit[1]"\t"region[1]"\t"region[2]}' $projPath/peakCalling/SEACR/${histName}_seacr_top0.01.peaks.stringent.bed >$projPath/peakCalling/SEACR/${histName}_seacr_top0.01.peaks.summitRegion.bed
computeMatrix reference-point -S $projPath/alignment/bigwig/${histName}_normalized_BPM_50bin.bw \
							   -p 24 \
							  --referencePoint center \
							  --afterRegionStartLength 3000 \
							  --beforeRegionStartLength 3000 \
							  -R $projPath/peakCalling/MACS2/macs2_${histName}_peak_q0.1_summits.bed \
							  --skipZeros \
							  -o $projPath/deeptools/${histName}_macs2_refPoint_center_data.gz

plotHeatmap -m $projPath/deeptools/${histName}_macs2_refPoint_center_data.gz \
            --missingDataColor 1 \
            --colorList 'white,#0066CC' \
            --heatmapHeight 12 \
			--sortUsing sum --startLabel "Peak Start" \
            --endLabel "Peak End" --xAxisLabel "" \
            --regionsLabel "Peaks" \
            --samplesLabel "${histName}" \
            -o $projPath/deeptools/${histName}_macs2_refPoint_center_heatmap.pdf
	}&
done		

##referencePoint center MACS2 narrowPeak
#0066CC 浅蓝色
ref_gtf="/storage/maxianjueLab/guoyifan/species_reference/Drosophila_melanogaster_BDGP6_32/Drosophila_melanogaster.BDGP6.32.109.gtf"
projPath="/storage/maxianjueLab/guoyifan/CUT_Tag"

for histName in E75_Flag1 E75_Flag2 E75_Flag4 ENY_N1 ENY_N2 ENY_N3 ENY_Y1 ENY_Y2 ENY_Y5 NY_N1 NY_N2 NY_N3 NY_Y1 NY_Y3 NY_Y4; do
{
awk '{print $1"\t"$2"\t"$3"\t"$4"\t+"}'  $projPath/peakCalling/MACS2/macs2_${histName}_peak_q0.1_peaks.narrowPeak > $projPath/peakCalling/MACS2/macs2_${histName}_peak_q0.1_narrowPeak.bed
computeMatrix reference-point -S $projPath/alignment/bigwig/${histName}_normalized_BPM_50bin.bw \
							   -p 24 \
							  --referencePoint center \
							  --afterRegionStartLength 3000 \
							  --beforeRegionStartLength 3000 \
							  -R $projPath/peakCalling/MACS2/macs2_${histName}_peak_q0.1_narrowPeak.bed \
							  --skipZeros \
							  -o $projPath/deeptools/${histName}_macs2_narrowPeak_center_data.gz

plotHeatmap -m $projPath/deeptools/${histName}_macs2_narrowPeak_center_data.gz \
            --missingDataColor 1 \
            --colorList 'white,#0066CC' \
            --heatmapHeight 12 \
			--sortUsing sum --startLabel "Peak Start" \
            --endLabel "Peak End" --xAxisLabel "" \
            --regionsLabel "Peaks" \
            --samplesLabel "${histName}" \
            -o $projPath/deeptools/${histName}_macs2_narrowPeak_center_heatmap.pdf
	}&
done		

# run compute matrix to collect the data needed for plotting
cd /storage/maxianjueLab/guoyifan/CUT_Tag/alignment/bigwig

##ENY_Y vs NY_Y   Yki Targets
computeMatrix scale-regions -S $projPath/alignment/bigwig/ENY_Y1_normalized_BPM_50bin.bw \
                             $projPath/alignment/bigwig/ENY_Y2_normalized_BPM_50bin.bw \
							 $projPath/alignment/bigwig/ENY_Y5_normalized_BPM_50bin.bw \
                             $projPath/alignment/bigwig/NY_Y1_normalized_BPM_50bin.bw \
							 $projPath/alignment/bigwig/NY_Y3_normalized_BPM_50bin.bw \
							 $projPath/alignment/bigwig/NY_Y4_normalized_BPM_50bin.bw \
                              -R $projPath/deeptools/targets/Yki_targets_bed/*.bed \
                              --beforeRegionStartLength 3000 \
                              --regionBodyLength 2000 \
                              --afterRegionStartLength 3000 \
                              --skipZeros -o $projPath/deeptools/targets/ENY_Y_vs_NY_Y_Yki_matrix.mat.gz

plotProfile -m $projPath/deeptools/targets/ENY_Y_vs_NY_Y_Yki_matrix.mat.gz \
     -out ENY_Y_vs_NY_Y_Yki_plotProfile.pdf \
     --plotType=fill \
     --perGroup \
     --colors red red red blue blue blue \
     --plotTitle "ENY_Y vs NY_Y Yki Targets"

##ENY_N vs NY_N Yki Targets
computeMatrix scale-regions -S $projPath/alignment/bigwig/ENY_N1_normalized_BPM_50bin.bw \
                             $projPath/alignment/bigwig/ENY_N2_normalized_BPM_50bin.bw \
							 $projPath/alignment/bigwig/ENY_N3_normalized_BPM_50bin.bw \
                             $projPath/alignment/bigwig/NY_N1_normalized_BPM_50bin.bw \
							 $projPath/alignment/bigwig/NY_N2_normalized_BPM_50bin.bw \
							 $projPath/alignment/bigwig/NY_N3_normalized_BPM_50bin.bw \
                              -R $projPath/deeptools/targets/Yki_targets_bed/*.bed \
                              --beforeRegionStartLength 3000 \
                              --regionBodyLength 2000 \
                              --afterRegionStartLength 3000 \
                              --skipZeros -o $projPath/deeptools/targets/ENY_N_vs_NY_N_Yki_matrix.mat.gz

plotProfile -m $projPath/deeptools/targets/ENY_N_vs_NY_N_Yki_matrix.mat.gz \
     -out ENY_N_vs_NY_N_Yki_plotProfile.pdf \
     --plotType=fill \
     --perGroup \
     --colors red red red blue blue blue \
     --plotTitle "ENY_N vs NY_N Yki Targets"

##ENY_Y vs NY_Y Notch Targets
computeMatrix scale-regions -S $projPath/alignment/bigwig/ENY_Y1_normalized_BPM_50bin.bw \
                             $projPath/alignment/bigwig/ENY_Y2_normalized_BPM_50bin.bw \
							 $projPath/alignment/bigwig/ENY_Y5_normalized_BPM_50bin.bw \
                             $projPath/alignment/bigwig/NY_Y1_normalized_BPM_50bin.bw \
							 $projPath/alignment/bigwig/NY_Y3_normalized_BPM_50bin.bw \
							 $projPath/alignment/bigwig/NY_Y4_normalized_BPM_50bin.bw \
                              -R $projPath/deeptools/targets/Notch_targets_bed/*.bed \
                              --beforeRegionStartLength 3000 \
                              --regionBodyLength 2000 \
                              --afterRegionStartLength 3000 \
                              --skipZeros -o $projPath/deeptools/targets/ENY_Y_vs_NY_Y_Notch_matrix.mat.gz

plotProfile -m $projPath/deeptools/targets/ENY_Y_vs_NY_Y_Notch_matrix.mat.gz \
     -out ENY_Y_vs_NY_Y_Notch_plotProfile.pdf \
     --plotType=fill \
     --perGroup \
     --colors red red red blue blue blue \
     --plotTitle "ENY_Y vs NY_Y Notch Targets"

##ENY_N vs NY_N Notch Targets
computeMatrix scale-regions -S $projPath/alignment/bigwig/ENY_N1_normalized_BPM_50bin.bw \
                             $projPath/alignment/bigwig/ENY_N2_normalized_BPM_50bin.bw \
							 $projPath/alignment/bigwig/ENY_N3_normalized_BPM_50bin.bw \
                             $projPath/alignment/bigwig/NY_N1_normalized_BPM_50bin.bw \
							 $projPath/alignment/bigwig/NY_N2_normalized_BPM_50bin.bw \
							 $projPath/alignment/bigwig/NY_N3_normalized_BPM_50bin.bw \
                              -R $projPath/deeptools/targets/Notch_targets_bed/*.bed \
                              --beforeRegionStartLength 3000 \
                              --regionBodyLength 2000 \
                              --afterRegionStartLength 3000 \
                              --skipZeros -o $projPath/deeptools/targets/ENY_N_vs_NY_N_Notch_matrix.mat.gz

plotProfile -m $projPath/deeptools/targets/ENY_N_vs_NY_N_Notch_matrix.mat.gz \
     -out ENY_N_vs_NY_N_Notch_plotProfile.pdf \
     --plotType=fill \
     --perGroup \
     --colors red red red blue blue blue \
     --plotTitle "ENY_N vs NY_N Notch Targets"


# make one image per BED file instead of per bigWig file
--perGroup \ 
# add color between the x axis and the lines
--plotType=fill \ 



















#6 生物学重复的处理

module load samtools/1.11
module load bedtools/2.30.0
module load bowtie/2.4.2
module load R/4.2.1
source ~/.bashrc
conda activate deeptools

#6.1 Merge bam files across biological replicates
cd /storage/maxianjueLab/guoyifan/CUT_Tag/alignment/bam

samtools merge E75_mapped.bam E75_Flag1_bowtie2.mapped.sorted.bam  E75_Flag2_bowtie2.mapped.sorted.bam  E75_Flag3_bowtie2.mapped.sorted.bam  E75_Flag4_bowtie2.mapped.sorted.bam &
samtools merge ENY_N_mapped.bam ENY_N10_bowtie2.mapped.sorted.bam  ENY_N2_bowtie2.mapped.sorted.bam  ENY_N6_bowtie2.mapped.sorted.bam  ENY_N8_bowtie2.mapped.sorted.bam &
samtools merge ENY_Y_mapped.bam ENY_Y2_bowtie2.mapped.sorted.bam  ENY_Y4_bowtie2.mapped.sorted.bam  ENY_Y6_bowtie2.mapped.sorted.bam  ENY_Y8_bowtie2.mapped.sorted.bam &
samtools merge NY_N_mapped.bam NY_N1_bowtie2.mapped.sorted.bam  NY_N5_bowtie2.mapped.sorted.bam  NY_N7_bowtie2.mapped.sorted.bam  NY_N9_bowtie2.mapped.sorted.bam &
samtools merge NY_Y_mapped.bam NY_Y1_bowtie2.mapped.sorted.bam  NY_Y3_bowtie2.mapped.sorted.bam  NY_Y5_bowtie2.mapped.sorted.bam  NY_Y7_bowtie2.mapped.sorted.bam &

projPath="/storage/maxianjueLab/guoyifan/CUT_Tag"

 for i in E75 NY_N NY_Y ENY_N ENY_Y;do
   {
     samtools index $projPath/alignment/bam/${i}_mapped.bam
	 bamCoverage -b $projPath/alignment/bam/${i}_mapped.bam -o $projPath/alignment/bigwig/${i}_raw.bw
	 bamCoverage -b $projPath/alignment/bam/${i}_mapped.bam -o $projPath/alignment/bigwig/${i}_normalized_BPM_50bin.bw --binSize 50 --normalizeUsing BPM	 
   }&
done

#暗绿色
for histName in E75 NY_N NY_Y ENY_N ENY_Y; do
{
computeMatrix reference-point -S $projPath/alignment/bigwig/${histName}_normalized_BPM_50bin.bw \
							   -p 24 \
							  --binSize 50 \
							  --regionBodyLength 5000 \
							  --referencePoint TSS \
							  --afterRegionStartLength 3000 \
							  --beforeRegionStartLength 3000 \
							  -R $ref_gtf \
							  --skipZeros \
							  -o $projPath/deeptools/${histName}_refPoint_TSS_TES_data.gz

plotHeatmap -m $projPath/deeptools/${histName}_refPoint_TSS_TES_data.gz \
            --missingDataColor 1 \
            --colorList 'white,#339933' \
            --heatmapHeight 12 \
			--sortUsing sum --startLabel "TSS" \
            --endLabel "TES" --xAxisLabel "" \
            --regionsLabel "Peaks" \
            --samplesLabel "${histName}" \
            -o $projPath/deeptools/${histName}_refPoint_TSS_TES_heatmap.pdf
    }&
done



#6.2 Overlapping peaks across biological replicates
bedtools intersect -h
projPath="/storage/maxianjueLab/guoyifan/CUT_Tag"
# 两个BED文件，第一步合并成一个文件，但来自两个文件的区域是分开的
cat $projPath/peakCalling/MACS2/macs2_E75_Flag1_peak_q0.1_narrowPeak.bed \
      $projPath/peakCalling/MACS2/macs2_E75_Flag2_peak_q0.1_narrowPeak.bed \
      $projPath/peakCalling/MACS2/macs2_E75_Flag4_peak_q0.1_narrowPeak.bed > $projPath/peakCalling/MACS2/E75_peak_q0.1_narrowPeak.bed

cat $projPath/peakCalling/MACS2/macs2_ENY_N1_peak_q0.1_narrowPeak.bed \
      $projPath/peakCalling/MACS2/macs2_ENY_N2_peak_q0.1_narrowPeak.bed \
      $projPath/peakCalling/MACS2/macs2_ENY_N3_peak_q0.1_narrowPeak.bed > $projPath/peakCalling/MACS2/ENY_N_peak_q0.1_narrowPeak.bed

cat $projPath/peakCalling/MACS2/macs2_NY_N1_peak_q0.1_narrowPeak.bed \
      $projPath/peakCalling/MACS2/macs2_NY_N2_peak_q0.1_narrowPeak.bed \
      $projPath/peakCalling/MACS2/macs2_NY_N3_peak_q0.1_narrowPeak.bed > $projPath/peakCalling/MACS2/NY_N_peak_q0.1_narrowPeak.bed	  

cat $projPath/peakCalling/MACS2/macs2_ENY_Y1_peak_q0.1_narrowPeak.bed \
      $projPath/peakCalling/MACS2/macs2_ENY_Y2_peak_q0.1_narrowPeak.bed \
      $projPath/peakCalling/MACS2/macs2_ENY_Y5_peak_q0.1_narrowPeak.bed > $projPath/peakCalling/MACS2/ENY_Y_peak_q0.1_narrowPeak.bed	  

cat $projPath/peakCalling/MACS2/macs2_NY_Y1_peak_q0.1_narrowPeak.bed \
      $projPath/peakCalling/MACS2/macs2_NY_Y3_peak_q0.1_narrowPeak.bed \
      $projPath/peakCalling/MACS2/macs2_NY_Y4_peak_q0.1_narrowPeak.bed > $projPath/peakCalling/MACS2/NY_Y_peak_q0.1_narrowPeak.bed

# 合并后的文件PEAK数是合并前两个文件的总和
# 执行合并之前，需要对染色体坐标进行排序，否则结果会有问题
# 排序后PEAK数不变
# 最后进行合并
# 这样即可得到我们需要的BED文件，PEAK数减少，但多于任何一个合并前的文件
cd $projPath/peakCalling/MACS2
 for i in E75 NY_N NY_Y ENY_N ENY_Y ;do
   {
sort -k1,1 -k2,2n ${i}_peak_q0.1_narrowPeak.bed > ${i}_peak_q0.1_narrowPeak_sort.bed
bedtools merge -d 10 -i ${i}_peak_q0.1_narrowPeak_sort.bed > ${i}_peak_q0.1_narrowPeak.bed
   }&
done

for histName in E75 NY_N NY_Y ENY_N ENY_Y; do
  {
macs2 callpeak -t ${projPath}/alignment/bam/${histName}_mapped.bam \
      -g dm -f BAMPE -n macs2_${histName}_peak_q0.1 --outdir $projPath/peakCalling/MACS2 -q 0.1 --keep-dup all 2>${projPath}/peakCalling/MACS2/macs2Peak_${histName}_summary.txt
    }&
done




##referencePoint center MACS2 narrowPeak
#0066CC 浅蓝色
ref_gtf="/storage/maxianjueLab/guoyifan/species_reference/Drosophila_melanogaster_BDGP6_32/Drosophila_melanogaster.BDGP6.32.109.gtf"
projPath="/storage/maxianjueLab/guoyifan/CUT_Tag"

for histName in E75 ENY_N ENY_Y NY_N NY_Y; do
{
computeMatrix reference-point -S $projPath/alignment/bigwig/${histName}_raw.bw \
							   -p 24 \
							  --referencePoint center \
							  --afterRegionStartLength 3000 \
							  --beforeRegionStartLength 3000 \
							  -R $projPath/peakCalling/MACS2/${histName}_peak_q0.1_narrowPeak.bed \
							  --skipZeros \
							  -o $projPath/deeptools/${histName}_macs2_narrowPeak_center_data.gz

plotHeatmap -m $projPath/deeptools/${histName}_macs2_narrowPeak_center_data.gz \
            --missingDataColor 1 \
            --colorList 'white,#0066CC' \
            --heatmapHeight 12 \
			--sortUsing sum --startLabel "Peak Start" \
            --endLabel "Peak End" --xAxisLabel "" \
            --regionsLabel "Peaks" \
            --samplesLabel "${histName}" \
            -o $projPath/deeptools/${histName}_macs2_narrowPeak_center_heatmap.pdf
	}&
done		




#6.3 Handling replicates using the Irreproducibility Discovery Rate (IDR) framework
conda create -n idr
conda install -c bioconda idr

source ~/.bashrc
conda activate idr
#Sort peak by -log10(p-value)
cd /storage/maxianjueLab/guoyifan/CUT_Tag/peakCalling/MACS2

ref_gtf="/storage/maxianjueLab/guoyifan/species_reference/Drosophila_melanogaster_BDGP6_32/Drosophila_melanogaster.BDGP6.32.109.gtf"
projPath="/storage/maxianjueLab/guoyifan/CUT_Tag"
for histName in E75_Flag1 E75_Flag2 E75_Flag4 ENY_N1 ENY_N2 ENY_N3 ENY_Y1 ENY_Y2 ENY_Y5 NY_N1 NY_N2 NY_N3 NY_Y1 NY_Y3 NY_Y4; do
{
sort -k8,8nr $projPath/peakCalling/MACS2/macs2_${histName}_peak_q0.1_peaks.narrowPeak > $projPath/peakCalling/MACS2/macs2_${histName}_IDR_peaks.narrowPeak
	}&
done	

mkdir -p 
for histName in E75_Flag1 E75_Flag2 E75_Flag4 ENY_N1 ENY_N2 ENY_N3 ENY_Y1 ENY_Y2 ENY_Y5 NY_N1 NY_N2 NY_N3 NY_Y1 NY_Y3 NY_Y4; do
cd /storage/maxianjueLab/guoyifan/CUT_Tag/peakCalling/IDR



idr --samples macs2_E75_Flag1_IDR_peaks.narrowPeak macs2_E75_Flag2_IDR_peaks.narrowPeak macs2_E75_Flag4_IDR_peaks.narrowPeak \
--input-file-type narrowPeak \
--rank p.value \
--output-file sample-idr \
--plot \
--log-output-file sample.idr.log














#7. Homer Annotation and Motif Enrichment
#find motif by homer
#https://mp.weixin.qq.com/s/mK13kmbblnBcu6aSSfsPTg
#https://mp.weixin.qq.com/s/tr2_zk1XvAyI1oHZPD6W-Q
#https://mp.weixin.qq.com/s/U7iAS75Mlg6aJFqGp1LfDg
wget http://homer.ucsd.edu/homer/configureHomer.pl

##使用configureHomer.pl配置Homer
conda create -n homer
conda install -c bioconda homer=4.11

perl configureHomer.pl -install  
perl /home/maxianjueLab/guoyifan/homer/configureHomer.pl -list
perl /home/maxianjueLab/guoyifan/homer/configureHomer.pl -install dm6
perl /home/maxianjueLab/guoyifan/miniconda3/envs/homer/share/homer/.//configureHomer.pl -install dm6

module load perl/5.34.0
module load samtools/1.11
module load bedtools/2.30.0
module load bowtie/2.4.2
module load R/4.2.1
projPath="/storage/maxianjueLab/guoyifan/CUT_Tag"
source ~/.bashrc
conda activate homer
#首先对bed文件进行处理,homer软件在读取bed文件时，需要提取对应的列作为输入文件。我们要对MACS找到的peaks记录文件，还需提取对应的列给HOMER作为输入文件，提取操作为：
#awk '{print $4"	"$1"	"$2"	"$3"	+"}' sample_peaks.bed >sample_homer.bed https://mp.weixin.qq.com/s/U7iAS75Mlg6aJFqGp1LfDg
#提取macs2输出文件中的sumit.bed为homer.bed
for histName in E75_Flag1 E75_Flag2 E75_Flag4 ENY_N1 ENY_N2 ENY_N3 ENY_Y1 ENY_Y2 ENY_Y5 NY_N1 NY_N2 NY_N3 NY_Y1 NY_Y3 NY_Y4; do
{
awk '{print $4"\t"$1"\t"$2"\t"$3"\t+"}' ${projPath}/peakCalling/MACS2/macs2_${histName}_peak_q0.1_peaks.narrowPeak > ${projPath}/peakCalling/homer/Positionsfile/macs2_${histName}_peak_q0.1_homer.tmp
  }&
done

for histName in E75 ENY_N ENY_Y NY_N NY_Y; do
{
awk '{print $4"\t"$1"\t"$2"\t"$3"\t+"}' ${projPath}/peakCalling/MACS2/macs2_${histName}_peak_q0.1_peaks.narrowPeak > ${projPath}/peakCalling/homer/Positionsfile/macs2_${histName}_peak_q0.1_homer.tmp
  }&
done

for histName in E75_Flag1 E75_Flag2 E75_Flag3 E75_Flag4 ENY_N10 ENY_N2 ENY_N6 ENY_N8 ENY_Y2 ENY_Y4 ENY_Y6 ENY_Y8 NY_N1 NY_N5 NY_N7 NY_N9 NY_Y1 NY_Y3 NY_Y5 NY_Y7; do
{
awk '{print $4"\t"$1"\t"$2"\t"$3"\t+"}' ${projPath}/peakCalling/MACS2/macs2_${histName}_peak_q0.1_peaks.narrowPeak > ${projPath}/peakCalling/homer/Positionsfile/macs2_${histName}_peak_q0.1_homer.tmp
  }&
done

#寻找富集Motifs
#findMotifsGenome.pl <Homer Peak/Positions file> <genome> <output directory> -size # [options]
ref="/storage/maxianjueLab/guoyifan/species_reference/Drosophila_melanogaster_BDGP6_32/Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa"
for histName in E75_Flag1 E75_Flag2 E75_Flag4 ENY_N1 ENY_N2 ENY_N3 ENY_Y1 ENY_Y2 ENY_Y5 NY_N1 NY_N2 NY_N3 NY_Y1 NY_Y3 NY_Y4; do
{
findMotifsGenome.pl ${projPath}/peakCalling/homer/Positionsfile/macs2_${histName}_peak_q0.1_homer.tmp ${ref} ${projPath}/peakCalling/homer/motifs/${histName} -len 8,10,12
  }&
done

for histName in E75 ENY_N ENY_Y NY_N NY_Y; do
{
findMotifsGenome.pl ${projPath}/peakCalling/homer/Positionsfile/macs2_${histName}_peak_q0.1_homer.tmp ${ref} ${projPath}/peakCalling/homer/motifs/${histName} -len 8,10,12,15,20
  }&
done
# 参数解释
#-输入文件：awk处理好的Homer Peak/Positions file
#-参考基因组：这里是hg19
#-输出文件：给一个路径和输出文件的名字
#-len：motif大小设置，默认8,10,12；越大需要的计算资源越多


#使用 annotatePeaks.pl 对peaks进行注释
ref_fa="/storage/maxianjueLab/guoyifan/species_reference/Drosophila_melanogaster_BDGP6_32/Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa"
ref_gtf="/storage/maxianjueLab/guoyifan/species_reference/Drosophila_melanogaster_BDGP6_32/Drosophila_melanogaster.BDGP6.32.109.gtf"
for histName in E75_Flag1 E75_Flag2 E75_Flag4 ENY_N1 ENY_N2 ENY_N3 ENY_Y1 ENY_Y2 ENY_Y5 NY_N1 NY_N2 NY_N3 NY_Y1 NY_Y3 NY_Y4; do
{
annotatePeaks.pl ${projPath}/peakCalling/homer/Positionsfile/macs2_${histName}_peak_q0.1_homer.tmp ${ref_fa} -gtf  ${ref_gtf} 1> ${projPath}/peakCalling/homer/Annotation/${histName}_peakAnn.xls 2> ${projPath}/peakCalling/homer/Annotation/${histName}_annLog.txt
  }&
done

for histName in E75 ENY_N ENY_Y NY_N NY_Y; do
{
annotatePeaks.pl ${projPath}/peakCalling/homer/Positionsfile/macs2_${histName}_peak_q0.1_homer.tmp ${ref_fa} -gtf  ${ref_gtf} 1> ${projPath}/peakCalling/homer/Annotation/${histName}_peakAnn.xls 2> ${projPath}/peakCalling/homer/Annotation/${histName}_annLog.txt
  }&
done

for histName in E75_Flag1 E75_Flag2 E75_Flag3 E75_Flag4 ENY_N10 ENY_N2 ENY_N6 ENY_N8 ENY_Y2 ENY_Y4 ENY_Y6 ENY_Y8 NY_N1 NY_N5 NY_N7 NY_N9 NY_Y1 NY_Y3 NY_Y5 NY_Y7; do
{
annotatePeaks.pl ${projPath}/peakCalling/homer/Positionsfile/macs2_${histName}_peak_q0.1_homer.tmp ${ref_fa} -gtf  ${ref_gtf} 1> ${projPath}/peakCalling/homer/Annotation/${histName}_peakAnn.xls 2> ${projPath}/peakCalling/homer/Annotation/${histName}_annLog.txt
  }&
done

#DiffBind差异peak分析
#bam文件排序https://mp.weixin.qq.com/s?__biz=MzUzNTM0NzM1Ng==&mid=2247483839&idx=1&sn=1619189704db5ff22442e9fad439e93f&chksm=fa87a84ccdf0215aff830c3081b236c1b1770bdcf1434626e1c10f96a00c44a21797cb543520&cur_album_id=2407970806904635394&scene=189#wechat_redirect
##linux代码


