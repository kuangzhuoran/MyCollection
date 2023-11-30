##PAM聚类
#参考链接如下
#https://mp.weixin.qq.com/s?__biz=Mzg3MzEwMzIwOA==&mid=2247484932&idx=1&sn=985a09909e18b1adc65cf48e99109d96&chksm=cee46c03f993e515a777d90cf5cdeac1a94701474853680b7dca5971fbe1e8432e3581a7f463&scene=21#wechat_redirect
#https://mp.weixin.qq.com/s?__biz=Mzg3MzEwMzIwOA==&mid=2247484997&idx=1&sn=55e844c1d4a695c9b44b5fb0ad5af3e6&chksm=cee46c42f993e55419b568b7c5b976f0ce35a578892bc9a3bb9a1c1ee482169ac52276267338&cur_album_id=1621136813132152835&scene=189#wechat_redirect
#https://mp.weixin.qq.com/s?__biz=Mzg3MzEwMzIwOA==&mid=2247485114&idx=1&sn=307cf1418063f89f84f9613dc8565f35&chksm=cee46cbdf993e5abd3286cbc7002d10b786d9584967247577edef491bf1d73a82916d0221173&cur_album_id=1621136813132152835&scene=189#wechat_redirect
#https://mp.weixin.qq.com/s?__biz=Mzg3MzEwMzIwOA==&mid=2247485150&idx=1&sn=15d5e5f2cb89ec99d553bd1f571272e5&chksm=cee46cd9f993e5cfeb742030efab86c25f74271f7676582ed2d384420703c7b059526d148249&cur_album_id=1621136813132152835&scene=189#wechat_redirect

library(vegan)
library(ggplot2)
library(fpc)
library(cluster)

data(dune)    #使用前务必查看这两个示例文件的数据格式
data(dune.env)
#dist = vegdist(a,method = "bray") 可以选择指定的距离
dune.dist <- vegdist(dune)

#pambest = pamk(dune.dist) #krange默认是2:10
pambest = pamk(dune.dist,krange = 2:19,criterion = 'ch')
#k-medoids聚类也就是PAM聚类
k=pambest$nc

pambest$crit   #这个就是不同K值的CHindex或者其他指标
#越高表示这个K值越好
#用这个画个图就行

#简单的画个图，不过这里的criterion只能是默认的asw
mypam = pam(dune.dist, k)
asw=numeric(nrow(dune))
for (i in 2:(length(asw)-1)) {
  asw[i]=pam(dune.dist, i)$silinfo$avg.width
}
k.best=which.max(asw)
plot(1:length(asw), asw, type="h", main="XD", lwd=2,
     xlab="k clusters", ylab="XD")
#简单的画个图，不过这里的criterion只能是默认的asw


##other
##k-means聚类
fit = cascadeKM(dune,2,19,iter=10,criterion="calinski") #请查看此函数的帮助文档
fit$results
calinski.best <- as.numeric(which.max(fit$results[2,]))
calinski.best

plot(fit, sortg = TRUE, grpmts.plot = TRUE)#这个画的右边那张图真tm反直觉

calinski<-as.data.frame(fit$results[2,])
calinski$cluster <- c(1:10)
ggplot(calinski,aes(x = calinski[,2], y = calinski[,1]))+geom_line()
