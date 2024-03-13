#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("ChIPseeker")
#https://mp.weixin.qq.com/s?__biz=MzI5NjUyNzkxMg==&mid=2247484887&idx=1&sn=899efd5476443f9c358b0b81214588c5&scene=21#wechat_redirect
library(ChIPseeker)
library(GenomicFeatures)
library(tidyverse)
#示例文件
files <- getSampleFiles()
peak2 <- readPeakFile(files[[4]])
head(peak2)

#基因注释
gff <- makeTxDbFromGFF('E:/workspace/spalax/0.data/golani.chr.v4.gff')
peak <- readPeakFile("conservation.CNE.bed")

peak <- readPeakFile('C:/Users/Administrator/Desktop/TmpFile/GY209-13.ATAC.LG08.bed')
peak <- readPeakFile('C:/Users/Administrator/Desktop/TmpFile/GY209-13.ATAC.other.bed')

gff <- makeTxDbFromGFF('C:/Users/Administrator/Desktop/TmpFile/Efontanierii.gff')
peak <- readPeakFile('C:/Users/Administrator/Desktop/TmpFile/ZH209-6.ATAC.bed')
#peak.bed，或者是别的任何bed
peak <- readPeakFile('E:/workspace/spalax/3DGenomics/TAD/galili.All.40k.TADboundaries.bed')
peak <- readPeakFile('before.chroder.judaei.largeINV.final.bed')

#第一列是染色体名，第二列第三列是起始/终止bp
#第四列是peak的名称，第五列是peak的值


peakAnno <- annotatePeak(peak, tssRegion = c(-3000,3000), TxDb = gff,
                         )

peakAnno.df <- as.data.frame(peakAnno)
write.table(peakAnno.df,sep = ',',"conservation.CNE.anno.csv")
peakAnno.gr <- as.GRanges(peakAnno)
head(peakAnno.gr, 3)

plotAnnoBar(peakAnno)#

vennpie(peakAnno)

plotAnnoPie(peakAnno)

plotDistToTSS(peakAnno)



#eg:TAD边界上在基因TSS上下游的分布
#但其实我想要的是基因在TAD边界上的分布
plotAvgProf2(peak, TxDb=gff, 
             upstream=3000, downstream=3000,
             xlab="Genomic Region (5'->3')", 
             ylab = "Read Count Frequency",
             conf = 0.95, resample = 1000)# 添加置信区间

#gff <- makeTxDbFromGFF('E:/workspace/spalax/0.data/galili.chr.v3.gff')
#peak.bed，或者是别的任何bed
#peak <- readPeakFile('E:/workspace/spalax/3DGenomics/TAD/galili.All.40k.TADboundaries.bed')
#第一列是染色体名，第二列第三列是起始/终止bp
#第四列是peak的名称，第五列是peak的值



#以下作废，还没有弄清楚

#  gff这个变量是可以其他的bed，比如TAD边界的bed
#那peak这个变量换成gene.bed
#那我们得到的就是gene在TAD边界上的分布
#换言之，是peak这个变量 在gff这个变量上的分布
gff <- makeTxDbFromGFF('E:/workspace/spalax/0.data/galili.chr.v3.gene.gff')
gff <- makeTxDbFromGFF('E:/workspace/spalax/0.data/galili.chr.v3.gff')

gff <- makeTxDbFromGFF('E:/workspace/SthIntereSting/ATAClike/galili.All.40k.TADboundaries.deal.gff')
peak <- readPeakFile('E:/workspace/SthIntereSting/ATAClike/galili.chr.v3.gene.bed')

plotAvgProf2(peak, TxDb=gff, 
             upstream=3000, downstream=3000,
             xlab="Genomic Region (5'->3')", 
             ylab = "Read Count Frequency",
             conf = 0.95, resample = 1000)# 添加置信区间

#启动子上下游3kb
promoter <- getPromoters(TxDb=gff, 
                         upstream=3000, downstream=3000)
#转为矩阵-peak在TSS上下游3kb的分布
tagMatrix <- getTagMatrix(peak, 
                          windows=promoter)
tagMatrix <- getTagMatrix(peak, 
                          windows=promoter,type = 'body')
#画成热图
tagHeatmap(tagMatrix, xlim=c(-3000, 3000), 
           color="red")


#多样本
go.gff <- makeTxDbFromGFF('E:/workspace/spalax/0.data/golani.chr.v4.gff')
ga.gff <- makeTxDbFromGFF('E:/workspace/spalax/0.data/galili.chr.v3.gff')
ca.gff <- makeTxDbFromGFF('E:/workspace/spalax/0.data/carmeli.chr.v4.gff')
ju.gff <- makeTxDbFromGFF('E:/workspace/spalax/0.data/judaei.chr.v3.gff')

golani = readPeakFile('golani.tmp1')
galili = readPeakFile('galili.tmp1')
carmeli= readPeakFile('carmeli.tmp1')
judaei = readPeakFile('judaei.tmp1')

write.table(goAnno.df,sep = ',',"golani.chrX.loop.anno.csv")
write.table(gaAnno.df,sep = ',',"galili.chrX.loop.anno.csv")
write.table(caAnno.df,sep = ',',"carmeli.chrX.loop.anno.csv")
write.table(juAnno.df,sep = ',',"judaei.chrX.loop.anno.csv")

golani = readPeakFile('E:/workspace/spalax/3DGenomics/TAD/golani.All.40k.TADboundaries.bed')
galili = readPeakFile('E:/workspace/spalax/3DGenomics/TAD/galili.All.40k.TADboundaries.bed')
carmeli= readPeakFile('E:/workspace/spalax/3DGenomics/TAD/carmeli.All.40k.TADboundaries.bed')
judaei = readPeakFile('E:/workspace/spalax/3DGenomics/TAD/judaei.All.40k.TADboundaries.bed')

go.gff <- makeTxDbFromGFF('~/TAD.boundary.gff')
golani <- readPeakFile('~/before.chroder.judaei.largeINV.final.sort.bed')
golani <- readPeakFile('~/INV.bp')

goAnno <- annotatePeak(golani, tssRegion = c(-3000, 3000), TxDb = go.gff)
goAnno.df <- as.data.frame(goAnno)
write.table(goAnno.df,sep = ',',"E:/workspace/spalax/3DGenomics/TAD/golani.All.40k.TADboundaries.anno.csv")

gaAnno <- annotatePeak(galili, tssRegion = c(-3000, 3000), TxDb = ga.gff)
gaAnno.df <- as.data.frame(gaAnno)
write.table(gaAnno.df,sep = ',',"E:/workspace/spalax/3DGenomics/TAD/galili.All.40k.TADboundaries.anno.csv")

caAnno <- annotatePeak(carmeli, tssRegion = c(-3000, 3000), TxDb = ca.gff)
caAnno.df <- as.data.frame(caAnno)
write.table(caAnno.df,sep = ',',"E:/workspace/spalax/3DGenomics/TAD/carmeli.All.40k.TADboundaries.anno.csv")

juAnno <- annotatePeak(judaei, tssRegion = c(-3000, 3000), TxDb = ju.gff)
juAnno.df <- as.data.frame(juAnno)
write.table(juAnno.df,sep = ',',"E:/workspace/spalax/3DGenomics/TAD/judaei.All.40k.TADboundaries.anno.csv")

peakAnnoList = list(S.golani=goAnno,S.galili=gaAnno,S.carmeli=caAnno,S.judaei=juAnno)
plotAnnoBar(peakAnnoList)
plotDistToTSS(peakAnnoList)
