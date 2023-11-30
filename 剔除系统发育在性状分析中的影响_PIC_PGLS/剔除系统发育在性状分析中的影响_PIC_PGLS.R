#系统发育独立差的实践 PIC
#示例文件同路径
#仔细浏览示例文件内容以及每一步代码带来的变化和效果
#来源https://mp.weixin.qq.com/s/GcI5vaBuPXvWaCS9Qn7kCA
#测试数据来自http://www.phytools.org/Cordoba2017/ex/3/PICs.html
obj<-read.csv("Centrarchidae.csv",row.names=1)
plot(obj$buccal.length, obj$gape.width)

fit.ols<-lm(gape.width~buccal.length,data=obj)
fit.ols
summary(fit.ols)
plot(obj$buccal.length, obj$gape.width)
abline(fit.ols,lwd=2,lty="dashed",col="red")

library(ape)
cent.tree<-read.tree("Centrarchidae.tre")

buccal.length<-setNames(obj[,"buccal.length"],
                        rownames(obj))
gape.width<-setNames(obj[,"gape.width"],rownames(obj))

#PIC处理
pic.bl<-pic(buccal.length,cent.tree)
pic.gw<-pic(gape.width,cent.tree)
fit.pic<-lm(pic.gw~pic.bl)

plot(pic.bl,pic.gw,xlab="PICs for buccal length",
     ylab="PICs for gape width",bg="grey",
     cex=1.4,pch=21)
abline(fit.pic,lwd=2,lty="dashed",col="red")
summary(fit.pic)


#如果使用随机生成的数据
library(phytools)
#随机生成一组数
set.seed(21) 
##生成随机树
tree<-rcoal(n=100)
plotTree(tree,ftype="off")
##模拟不具有相关性的布朗运动
x<-fastBM(tree)
y<-fastBM(tree)
par(mar=c(5.1,4.1,2.1,2.1))
plot(x,y,cex=1.4,pch=21,bg="grey")
fit<-lm(y~x)
abline(fit,lwd=2,lty="dashed",col="red")
fit
summary(fit)
###p-value: 0.003809
anova(fit)
###注意：以上随机生成的数据并不具有进化上的相关性
###但是拟合的结果显示P值显著
###显然假设检验出现了假阳性的错误


phylomorphospace(tree,cbind(x,y),label="off",node.size=c(0,0))
points(x,y,pch=21,bg="grey",cex=1.4)
abline(fit,lwd=2,lty="dashed",col="red")
##进行PIC的处理
ix<-pic(x,tree)
iy<-pic(y,tree)
fit<-lm(iy~ix-1) ##不明白为什么去掉截距项
fit
summary(fit)
###p-value: 0.6567
anova(fit)
###经过PIC处理之后的数据，不具有显著的相关性





#系统发育广义最小二乘法PGLS
#和上面的系统发育独立差PIC一样，PGLS也是为了在性状分析中，减少由系统发育关系带来的假阳性错误。
#测试数据来自http://www.phytools.org/Cordoba2017/ex/4/PGLS.html
#包含为研究鸟类鸣声进化，采集的33个性状数据；以及物种进化树文件。
library(ape)
library(nlme)
library(geiger)

#1：数据准备
##读取性状数据
datos<-read.csv("Barbetdata.csv",header=TRUE,row.names=1)
datos
##读取进化树
arbol<-read.nexus("BarbetTree.nex.txt")
arbol   
##读取为普通phylo格式，有42个末端节点，也就是样本
##检查性状数据与进化树中的物种是否匹配
obj<-name.check(arbol,datos)
obj
##进化树中存在9个性状数据中没有的样本
##去除不匹配的末端节点

arbol.cortado<-drop.tip(arbol, obj$tree_not_data)
name.check(arbol.cortado,datos)
##如果是进化树和性状各有一些不匹配的情况
##可以直接使用picante程序包中match.phylo.data()
#combined <- match.phylo.data(phy,trait)
##用combined中的数据，替换原始数据即可。
#phy <- combined$phy
#trait <- combined$data


##普通gls作为对照
modelo0 <- gls(Lnote~Lnalt, data=datos)
summary(modelo0)
anova(modelo0)
#Denom. DF: 31
#numDF  F-value p-value
#(Intercept)     1 34.15493  <.0001
#Lnalt           1  2.00933  0.1663   #一般关注这个P值
##Intercept的P值一般不关注，为满足误差项N(0,1)设置
##因而关注变量X的P-value即可，变量P值不显

##pgls
modelo1<-gls(Lnote~Lnalt, data=datos, correlation=bm)
summary(modelo1)
anova(modelo1)
?gls
