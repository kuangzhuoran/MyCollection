###距离矩阵的相关性(linkET)


# https://github.com/Hy4m/linkET  参考链接
#链接到 github 中安装 linkET 包。
#install.packages('devtools')
#devtools::install_github('Hy4m/linkET')

library(linkET)
library(ggplot2)
library(dplyr)

#https://mp.weixin.qq.com/s/bqEYio4RTbGwtO0n5iu1vw


###Mantel test
#和上文区别在于，Mantel tests是直接对两个距离矩阵进行相关性分析
#而LinkET则是从两种组成表(物种分度表和环境因子)开始分析
#当然这是指输入文件的不同而已

library(vegan)

#输入数据
#数据格式是              ---------(1)
#Sample OTU1 OTU2 OTU3   
#  A      1  0.8  0.9
#  B      2  0.7  0.6 
MT1 <- read.csv("MT1.csv", header = T, row.names = 1, check.names = F)
Gene1 <- read.csv("Gene1.csv", header = T, row.names = 1, check.names = F)
#计算输入数据的距离矩阵
MT1_DM <- vegdist(MT1, method = "bray")
Gene1_DM <- vegdist(Gene1, method = "euclidean")
#这里，对于多维数据的距离矩阵计算，使用bray-curtis算法；对于单变量（一维）数据，使用euclidean算法
#执行mantel test检验
set.seed(520) #设置随机种子，保证相同数据每次执行结果都一样
mantel(MT1_DM, Gene1_DM, method="spearman", permutations=9999) #计算MT1类群与Gene1的相关性
##这里，置换检验9999次，可更改以比较结果差异


#如果直接输入距离矩阵则是这样的数据格式
#     B1  B2 C1  C2
# B1  1   2  1   3
# B2  2   4  3   7 
soilgene = read.table("sg.tsv",sep = " ", header = T, row.names = 1)
gutgene =  read.table("gg.tsv",sep = " ", header = T, row.names = 1)
#bray的获取如下
gut_species_rate_raw = read.csv("E:/workspace/metagenome/gut/kraken2_reads/brackenReport0/raw/species.tmp.csv", header = T, row.names = 1)
#这个文件即 （1） 的转置版
soil_species_rate_raw = read.csv("E:/workspace/metagenome/soil/kraken2/raw/species.tmp.csv",header = T, row.names = 1)

t_gut_species_rate_raw = t(gut_species_rate_raw)
gutbray = vegdist(t_gut_species_rate_raw,method = "bray")
gutbray <- as.matrix(gutbray)

t_soil  = t(soil_species_rate_raw)
soilbray = vegdist(t_soil,method = "bray")
soilbray = as.matrix(soilbray)

# mantel(gutbray, soilbray, method="spearman", permutations=9999)
mantel(soilgene, gutgene, method="spearman", permutations=9999)
mantel(soilbray, gutgene, method="spearman", permutations=9999)

library(ape)

tr = bionj(bray)
plot(tr)

tr.nwk = write.tree(tr)

outdist = dist(bray)
out.hclust = hclust(outdist)
hclusttree = as.phylo(out.hclust)

hclust.nwk = write.tree(hclusttree)

mantel(host, bray, method="spearman", permutations=9999)


