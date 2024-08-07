library(nlme)
library(lme4)
library(vegan)
library(tidyverse)

setwd("D:/KZR/note_2021_5_19/线性混合模型")

#bp_data <- read.csv("ChrVNE.2Groups.csv",header = T)
bp_data <- read.csv("ChrVNE.csv",header = T)
# 检查是否有NA值
anyNA(bp_data)
# 查看前几行数据
head(bp_data)

mod1 <- lme(VNE ~ group + length, random = ~1 | id, data = bp_data)
mod1 <- lme(VNE ~ group, random = ~1 | species, data = bp_data)

mod1 <- lmer(VNE ~ group*length + (1|species), data = bp_data)

#此时固定效应是group(即是否治疗)、time(治疗的时间)、以及二者的交互
#随机效应为id,即不同的重复/观测值
summary(mod1)
anova(mod1)

# 增加交互效应
mod2 <- update(mod1, .~. + group:length)
summary(mod2)

anova(mod2)
anova(mod1,mod2)

# 增加随机斜率
mod3 <- update(mod1, random = ~1 + species | id)
summary(mod3)
anova(mod3)
anova(mod1,mod3)
anova(mod2,mod3)



##
library(vegan)
t = aov(bp_data$VNE~bp_data$length*bp_data$group)
t
summary(t)

##
library(GLMMadaptive)
setwd("D:/KZR/note_2021_5_19/线性混合模型")
bp_data <- read.csv("ChrVNE.csv",header = T)

fm <- mixed_model(fixed = VNE ~ group + length, random = ~ 1 | id, data = bp_data,
                  family = poisson())

