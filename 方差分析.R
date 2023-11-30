library(vegan)

#ADONIS，被称为多元方差分析，亦称为PERMANOVA (Permutational multivariate analysis of variance，置换多元方差分析)或非参数多因素方差分析(nonparametric MANOVA)，
#是一种基于样品距离（默认为Bray-Curtis或者其余距离如Euclidean)的非参数多元方差分析，是MANOVA的等同形式。
#它通过线性模型分析不同分组因素或环境因子(如土壤理化性质等)对样品差异的解释度，并使用置换检验进行显著性分析

data(dune)
data(dune.env)
dune.div <- adonis2(dune ~ A1, data = dune.env, permutations = 999, method="bray")
#看A1这个环境因素的影响是否显著

a = read.csv("tmp.csv", header = T, row.names = 1)
a.env = read.csv("tmp.env.csv", header = T, row.names = 1)
a.div <- adonis2(a ~ spp, data = a.env, permutations = 999, method="bray")

#如果这里不加by = "margin"参数，那么把A1和Management互换顺序，结果也不同，表明不同环境因素之间存在相关性
#为了排除这种相关性，使用这个参数：
adonis2(dune ~ A1 + Management, data = dune.env, permutations = 999, method="bray", by = "margin")
adonis2(dune ~ Management + A1, data = dune.env, permutations = 999, method="bray", by = "margin")


CaZyLv2 = read.csv("E:/workspace/metagenome/gut/MAGs_Bin/Per/all.coverm.csv", header = T, row.names = 1)
tCaZyLv2 = t(CaZyLv2)
env = read.csv("E:/workspace/metagenome/gut/MAGs_Bin/Per/1.csv", header = T)

div = adonis2(tCaZyLv2 ~ group, data = env, permutations = 1000, method = "bray") #??????(group)????????

#检验数据离散程度
bray = vegdist(tCaZyLv2,method = "bray") #???þ???????
dispersion = betadisper(bray, group=env$group)
permutest(dispersion) #??ɢ?ȼ???
#如果结果不显著，表明不同分组的离散度即方差没有显著差异
#即我们用方差分析检验出的差异是因为均值不同造成的而不是离散度 (每组数据在空间的中心点)


##一元单方差分析 #ANOVA
#多次的一元方差分析 VS 一次的多元方差分析

#元即因变量，多次的一元方差分析要求因变量之间互相独立

#https://zhuanlan.zhihu.com/p/296726829
data(iris)   #注意此数据和前文adonis示例数据的异同


aov1 <- aov(Sepal.Length~Species, iris) 
#其中Species将鸢（yuan1）尾花分为三种，其它4个指标分别表示鸢尾花的花的特征。比如Sepal.Length表示花瓣的长度。
#这里Species成为分组变量，也就是所谓的“因素”（相当于回归分析中的x）；
#花瓣的长度为结果变量（响应变量，相当于回归分析中的y），是比较组间均值差异的变量指标。

aov1
summary(aov1)
#结果解释：p值远小于0.001，说明鸢尾花的品种这个因素，对鸢尾花花瓣的长度，有显著的影响；
#或通俗点说，不同种类的鸢尾花的花瓣长度（均值）显著不同（至少有两种显著不同）

#为了进一步得到更详细的两两比较，我们可以使用R中的TukeyHSD()函数实现
tukey <- TukeyHSD(aov1)
tukey = as.data.frame(tukey$Species)

tukey$pair = rownames(tukey)
##多元方差分析参考以下链接
#https://mp.weixin.qq.com/s?__biz=Mzg3MzEwMzIwOA==&mid=2247484426&idx=1&sn=de99f5ef6dd069e16ea8898b5973fb50&chksm=cee46e0df993e71b215263e859afd481386654b74af505105ed12f1b27eba84cc1c54a4bf03c&scene=178&cur_album_id=1870267717979537411#rd
aov2 <- aov(Sepal.Length~Petal.Length*Sepal.Width, iris) 
aov2
summary(aov2)
