#组间差异分析：Anosim

#Anosim分析（Analysis of similarities）是一种基于置换检验和秩和检验的非参数检验方法
#用来检验组间的差异是否显著大于组内差异，从而判断分组是否有意义
#https://www.jianshu.com/p/dfa689f7cafd 这个链接是直接用16S跑的示例

#ANOSIM，相似性分析是一种非参数检验，用于检验高纬度数据间的相似性
#比较组间和组内差异的大小，从而判断分组是否有意义
#其可以用于检验两组的组间和组内差异，也可以用于多组
#https://www.jianshu.com/p/dfa689f7cafd

library(vegan)
library(ggplot2)
data(dune)    #使用前务必查看这两个示例文件的数据格式
data(dune.env)

# 首先需要将数据转成dist格式
#http://www.360doc.com/content/21/1128/23/73584707_1006323091.shtml
#dist = vegdist(a,method = "bray") 可以选择指定的距离
dune.dist <- vegdist(dune)
# 然后再利用anosim函数进行计算
dune.ano <- anosim(dune.dist, dune.env$Management, permutations = 999)
summary(dune.ano)

#ANOSIM statistic R为统计量，分布衡量的为零模型的分布
#Upper quantiles of permutations就是通过999次置换得到统计量的分位数
#如果R>0，说明组内距离小于组间距离，即分组是有效的

#可以进一步提取距离的秩的分析结果，并进行绘制：
dune.ano$dis.rank
dune.ano$class.vec
plot(dune.ano)


#从结果来看，可以看出NM分组效果较差，但总体来说分组是有效的
#简单做个降维，也能发现NM的点更分散
CaZybray = as.matrix(dune.dist)
pcoa = cmdscale(CaZybray,k=10,eig = T)
poi = pcoa$points
eigval = pcoa$eig
pcoa_eig = (pcoa$eig)[]/sum(pcoa$eig)  
poi = as.data.frame(poi)

poi2 = cbind(poi, dune.env$Management)

pcoa_plot = ggplot(data = poi2,aes(x=V1,y=V2))+  theme_bw() + 
  theme(panel.grid = element_blank()) +
  geom_point(aes(color=dune.env$Management,shape=dune.env$Management),size=4)



