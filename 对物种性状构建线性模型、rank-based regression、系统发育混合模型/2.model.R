####
##0.这里是一个LMM的示例代码 - 请忽略
####
library(nlme)
#library(lme4)
library(vegan)
library(tidyverse)
library(lmerTest)

setwd("E:/workspace/vision/17.trait.model")
bp_data <- read.csv("test.VNE.csv",header = T)
head(bp_data)
mod7 <- lme(VNE ~ group + length + species, random = ~1 | id, data = bp_data)
summary(mod7)

####
#1. 有无系统发育信号
####
library(ape)
library(phytools)
library(dplyr)
library(tibble)

## 读入数据
dat2 <- read.csv("test.VA.AL.celldensity.csv", header = TRUE)
tree2 <- read.tree("subset_tree.nwk")

## 按树 tip 顺序重排数据
dat2 <- dat2[match(tree2$tip.label, dat2$species), ]

## 再确认一次顺序完全一致
stopifnot(all(tree2$tip.label == dat2$species))

## 可选：对 RGC peak density 再做一个 log10 版本
#dat2$log10_RGC.peak.density <- log10(dat2$RGC.peak.density)

## 定义一个函数，自动跑 K 和 lambda
get_phylosig <- function(trait, tree, dat){
  x <- dat[[trait]]
  names(x) <- dat$species
  
  ## 去掉 NA
  keep <- !is.na(x)
  x <- x[keep]
  tree_sub <- keep.tip(tree, names(x))
  
  ## 按树顺序重排
  x <- x[tree_sub$tip.label]
  
  ## K
  resK <- phylosig(tree_sub, x, method = "K", test = TRUE)
  
  ## lambda
  resL <- phylosig(tree_sub, x, method = "lambda", test = TRUE)
  
  tibble(
    trait = trait,
    n_species = length(x),
    K = unname(resK$K),
    P_K = unname(resK$P),
    lambda = unname(resL$lambda),
    P_lambda = unname(resL$P)
  )
}

## 依次检验
traits_to_test <- c("maxVA", "AL", "RGC.peak.density", "log10_RGC.peak.density")

signal_results <- bind_rows(lapply(traits_to_test, get_phylosig, tree = tree2, dat = dat2))
#
#运行到这里就可以了
print(signal_results)
#
## 如果你还想看每个性状的详细输出
res_maxVA_K <- phylosig(tree2, setNames(dat2$maxVA, dat2$species), method = "K", test = TRUE)
res_maxVA_L <- phylosig(tree2, setNames(dat2$maxVA, dat2$species), method = "lambda", test = TRUE)

res_AL_K <- phylosig(tree2, setNames(dat2$AL, dat2$species), method = "K", test = TRUE)
res_AL_L <- phylosig(tree2, setNames(dat2$AL, dat2$species), method = "lambda", test = TRUE)

res_RGC_K <- phylosig(tree2, setNames(dat2$RGC.peak.density, dat2$species), method = "K", test = TRUE)
res_RGC_L <- phylosig(tree2, setNames(dat2$RGC.peak.density, dat2$species), method = "lambda", test = TRUE)

#res_logRGC_K <- phylosig(tree2, setNames(dat2$log10_RGC.peak.density, dat2$species), method = "K", test = TRUE)
#res_logRGC_L <- phylosig(tree2, setNames(dat2$log10_RGC.peak.density, dat2$species), method = "lambda", test = TRUE)

res_maxVA_K
res_maxVA_L
res_AL_K
res_AL_L
res_RGC_K
res_RGC_L
#res_logRGC_K
#res_logRGC_L

####
#3. PIC + lm -  如果残差满足正态分布
####
library(ape)
## 读入数据
dat2 <- read.csv("test.VA.AL.celldensity.csv", header = TRUE)
tree2 <- read.tree("subset_tree.nwk")
dat2 <- dat2[match(tree2$tip.label, dat2$species), ]
#dat2$log10_RGC <- log10(dat2$RGC.peak.density)

pic_maxVA <- pic(setNames(dat2$maxVA, dat2$species), tree2)
pic_AL <- pic(setNames(dat2$AL, dat2$species), tree2)
pic_RGC <- pic(setNames(dat2$RGC.peak.density, dat2$species), tree2)
#pic_log10_RGC <- pic(setNames(dat2$log10_RGC, dat2$species), tree2)

###一些数据诊断
par(mfrow = c(2, 2))
plot(mod_full)

mod_full <- lm(pic_maxVA ~ pic_AL + pic_RGC - 1)
res <- resid(mod_full)

qqnorm(res)
qqline(res)

hist(res, breaks = 20)
shapiro.test(res)

mod_full <- lm(pic_maxVA ~ pic_AL + pic_RGC - 1)
summary(mod_full)
#第一，在完整模型里，RGC 的统计信号更强。
#pic_AL   t = 4.461
#pic_RGC  t = 5.952
#在同一个 PIC 多元回归里，控制另一个变量后，
#绝对值更大的 t 值通常说明这个变量的独立解释力更强。这里是 RGC > AL。

#####

mod_AL_only <- lm(pic_maxVA ~ pic_AL - 1)
summary(mod_AL_only)

mod_RGC_only <- lm(pic_maxVA ~ pic_RGC - 1)
summary(mod_RGC_only)

#第二，看删减模型时，去掉 RGC 带来的模型损失更大
anova(mod_AL_only, mod_full)
#表示在只有 AL 的基础上，再加入 RGC，模型改进了：
#Sum of Sq = 12490
#F = 35.423
#P = 8.83e-08
anova(mod_RGC_only, mod_full)
#Sum of Sq = 7015.6
#F = 19.896
#P = 2.95e-05

#RGC 带来的额外解释量比 AL 更大
#第三，看单变量模型时，RGC 单独就能解释更多变异。
#mod_AL_only：
#R-squared = 0.1254
#mod_RGC_only：
#R-squared = 0.2518


#偏相关系数
tt <- coef(summary(mod_full))[,"t value"]
df <- df.residual(mod_full)

pr2_from_t <- tt^2 / (tt^2 + df)
pr2_from_t


####
#4. PIC + rank based regression -  如果残差不满足正态分布
####
library(ape)
library(Rfit)
dat2 <- read.csv("test.VA.AL.celldensity.csv", header = TRUE)
tree2 <- read.tree("subset_tree.nwk")
dat2 <- dat2[match(tree2$tip.label, dat2$species), ]
#dat2$log10_RGC <- log10(dat2$RGC.peak.density)

pic_maxVA <- pic(setNames(dat2$maxVA, dat2$species), tree2)
pic_AL <- pic(setNames(dat2$AL, dat2$species), tree2)
pic_RGC <- pic(setNames(dat2$RGC.peak.density, dat2$species), tree2)
#pic_log10_RGC <- pic(setNames(dat2$log10_RGC, dat2$species), tree2)


mod_rfit_full <- rfit(pic_maxVA ~ pic_AL + pic_RGC - 1)
summary(mod_rfit_full)
#在完整模型里，看 AL 和 RGC 的系数、检验统计量和 P 值
#如果二者都显著，说明在 rank based regression 下也成立
#谁的统计证据更强，可以先看 full model 里的检验统计量和 P 值

#Estimate
#这是 rank based regression 的回归系数估计。
#pic_AL = 0.66126
#pic_RGC = 9.8805e-05
#解释方向和普通回归一样：
#pic_AL 的系数为正，说明在控制 RGC 后，AL 的独立对比越大，maxVA 的独立对比也倾向越大
#pic_RGC 的系数也为正，说明在控制 AL 后，RGC 的独立对比越大，maxVA 的独立对比也倾向越大#
#也就是说，两者和 maxVA 都是正相关。
#但这里仍然不能直接比较两个 Estimate 的绝对值大小，因为 AL 和 RGC 的量纲完全不同。
#系数数值本身不能直接回答“谁影响更大”。

#Std. Error#
#这是 rank based regression 下对应系数的标准误。
#AL 的标准误 0.0682
#RGC 的标准误 9.44e-06
#标准误越小，表示该系数估计越稳定。

#t.value
#这是 rank based regression 的 t 型检验统计量。
#AL: 9.6937
#RGC: 10.4643
#在同一个完整模型里，绝对值更大的 t.value 通常表示独立统计证据更强。所以按这个标准，还是：
#RGC 略强于 AL。

#p.value
#两个变量都极显著：
#AL: 1.243e-14
#RGC: 4.945e-16
#这说明在 rank based regression 下，
#即使不依赖普通 lm 的正态残差假设，AL 和 RGC 仍然都显著。

#Robust R-squared 可以粗略理解为：
#在 rank based regression 框架下，这个模型对响应变量变异的解释强度

#Multiple R-squared (Robust): 0.6255977 
#这里你可以把它理解成：

#rank based regression 版的“整体模型检验”
#零假设大致就是：
#AL 和 RGC 对 pic_maxVA 都没有解释作用
#
#检验统计量很大
#
#p-value: 0
#这里的 0 不是数学上真等于 0，而是说明 P 值小到打印成 0 了。也就是说：
#整体模型极显著。

#####

mod_rfit_AL_only <- rfit(pic_maxVA ~ pic_AL - 1)
summary(mod_rfit_AL_only)
#Coefficients:
#  Estimate Std. Error t.value   p.value    
#pic_AL 0.535575   0.093301  5.7403 2.088e-07 ***

#单独只放 AL 时，AL 仍然显著正相关
#t.value = 5.7403
#P = 2.088e-07
#即使不考虑 RGC，AL 单独对 maxVA 也是有显著解释力的。
#Robust R-squared = 0.4154，说明单独 AL 的 rank based 模型也有一定解释力。

mod_rfit_RGC_only <- rfit(pic_maxVA ~ pic_RGC - 1)
summary(mod_rfit_RGC_only)

#Coefficients:
#  Estimate Std. Error t.value   p.value    
#pic_RGC 9.8188e-05 1.4085e-05  6.9713 1.248e-09 ***

#单独只放 RGC 时，RGC 也显著正相关
#t.value = 6.9713
#P = 1.248e-09
#而且从单变量模型看，RGC 的统计证据比 AL 更强：

#RGC only: t = 6.97
#AL only: t = 5.74
#所以单看一对一关系，RGC 也比 AL 更强。


####
#4. PGLS / phylogenetic linear regression + 其实还是需要残差满足正态的
####
library(ape)
library(phylolm)
## 读入数据
dat2 <- read.csv("test.VA.AL.celldensity.csv", header = TRUE)
tree2 <- read.tree("subset_tree.nwk")
dat2 <- dat2[match(tree2$tip.label, dat2$species), ]
stopifnot(all(tree2$tip.label == dat2$species))

rownames(dat2) <- dat2$species

mod_phy_full <- phylolm(maxVA ~ AL + RGC.peak.density,data = dat2,phy = tree2,
                        model = "lambda")
summary(mod_phy_full)

#lambda = 1
#说明系统发育相关性非常强，达到lambda模型上限

#AL 的系数是正的，而且显著：
#Estimate = 0.6150
#P = 2.95e-05
#说明在控制 RGC 后，AL 越大，maxVA 越高。

#RGC.peak.density 的系数也是正的，而且更显著：
#Estimate = 1.1345e-04
#P = 8.83e-08
#说明在控制 AL 后，RGC 越大，maxVA 越高。

#模型解释度
#R-squared = 0.4138
#AL 和 RGC 联合起来解释了大约 41.4% 的 maxVA 变异

#####
#删减模型
#####
mod_phy_AL_only <- phylolm(maxVA ~ AL,data = dat2,phy = tree2,model = "lambda")
summary(mod_phy_AL_only)

mod_phy_RGC_only <- phylolm(maxVA ~ RGC.peak.density,data = dat2,phy = tree2,model = "lambda")
summary(mod_phy_RGC_only)

#####
#比较删掉谁对模型影响更大
#####

## 1. 看 AIC
AIC(mod_phy_full)
AIC(mod_phy_AL_only)
AIC(mod_phy_RGC_only)
#看 AIC
#full: 471.96
#AL only: 499.96
#RGC only: 488.26
#所以：
#去掉 RGC 后，AIC 增加了 499.96 - 471.96 = 28.01
#去掉 AL 后，AIC 增加了 488.26 - 471.96 = 16.30
#去掉 RGC 让模型变差更多。
#而且两个 ΔAIC 都大于 10，说明这两个单变量模型都明显比 full model 差。

## 2. 看 logLik
logLik(mod_phy_full)
logLik(mod_phy_AL_only)
logLik(mod_phy_RGC_only)
## 3. 似然比检验
## full vs 去掉 RGC
LR_RGC <- 2 * (as.numeric(logLik(mod_phy_full)) - as.numeric(logLik(mod_phy_AL_only)))
P_RGC <- pchisq(LR_RGC, df = 1, lower.tail = FALSE)

## full vs 去掉 AL
LR_AL <- 2 * (as.numeric(logLik(mod_phy_full)) - as.numeric(logLik(mod_phy_RGC_only)))
P_AL <- pchisq(LR_AL, df = 1, lower.tail = FALSE)

LR_RGC
P_RGC

LR_AL
P_AL
#去掉 RGC，LR 更大，P 更小

####
#dat2$AL_z <- scale(dat2$AL)
#dat2$RGC_z <- scale(dat2$RGC.peak.density)
#mod_phy_full_z <- phylolm(maxVA ~ AL_z + RGC_z,data = dat2,phy = tree2,model = "lambda")
#summary(mod_phy_full_z)

###
###
#之前的model都选的lambda
#那如何考量要不要选其他model呢？
#跑几个模型，比AIC即可
mod_BM <- phylolm(maxVA ~ AL + RGC.peak.density,
                  data = dat2, phy = tree2, model = "BM")
mod_lambda <- phylolm(maxVA ~ AL + RGC.peak.density,
                      data = dat2, phy = tree2, model = "lambda")
mod_kappa <- phylolm(maxVA ~ AL + RGC.peak.density,
                     data = dat2, phy = tree2, model = "kappa")
mod_delta <- phylolm(maxVA ~ AL + RGC.peak.density,
                     data = dat2, phy = tree2, model = "delta")

model_comp <- data.frame(
  model = c("BM", "lambda", "kappa", "delta"),
  AIC = c(AIC(mod_BM),
          AIC(mod_lambda),
          AIC(mod_kappa),
          AIC(mod_delta)),
  logLik = c(as.numeric(logLik(mod_BM)),
             as.numeric(logLik(mod_lambda)),
             as.numeric(logLik(mod_kappa)),
             as.numeric(logLik(mod_delta)))
)

model_comp <- model_comp[order(model_comp$AIC), ]
model_comp$deltaAIC <- model_comp$AIC - min(model_comp$AIC)
model_comp
#不同模型之间的AIC差值很小，小于2
#所以选lambda也没关系
#

###
###


####
#5. PGLS / phylogenetic linear regression + bootstrap
#我暂时没太懂这个和残差正态之间的关系
####
library(ape)
library(phylolm)

## 读入数据
dat2 <- read.csv("test.VA.AL.celldensity.csv", header = TRUE)
tree2 <- read.tree("subset_tree.nwk")
dat2 <- dat2[match(tree2$tip.label, dat2$species), ]
stopifnot(all(tree2$tip.label == dat2$species))
rownames(dat2) <- dat2$species

## bootstrap 版主模型
mod_phy_full_boot <- phylolm(maxVA ~ AL + RGC.peak.density,
                             data = dat2,
                             phy = tree2,
                             model = "lambda",
                             boot = 1000,
                             full.matrix = TRUE)

summary(mod_phy_full_boot)

## bootstrap 结果
mod_phy_full_boot$bootmean
mod_phy_full_boot$bootsd
mod_phy_full_boot$bootconfint95

## 如果想看全部 bootstrap 抽样结果
#head(mod_phy_full_boot$bootstrap)

#####
#删减模型
#####
mod_phy_AL_only_boot <- phylolm(maxVA ~ AL,data = dat2,phy = tree2,model = "lambda",boot = 1000)
summary(mod_phy_AL_only_boot)

mod_phy_RGC_only_boot <- phylolm(maxVA ~ RGC.peak.density,data = dat2,phy = tree2,model = "lambda",boot = 1000)
summary(mod_phy_RGC_only_boot)

## 1. 看 AIC
AIC(mod_phy_full_boot)
AIC(mod_phy_AL_only_boot)
AIC(mod_phy_RGC_only_boot)
#看 AIC
#full: 471.96
#AL only: 499.96
#RGC only: 488.26
#所以：
#去掉 RGC 后，AIC 增加了 499.96 - 471.96 = 28.01
#去掉 AL 后，AIC 增加了 488.26 - 471.96 = 16.30
#去掉 RGC 让模型变差更多。
#而且两个 ΔAIC 都大于 10，说明这两个单变量模型都明显比 full model 差。
