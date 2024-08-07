library(nlme)
#library(lme4)
library(vegan)
library(tidyverse)
library(lmerTest)

setwd("D:/KZR/note_2021_5_19/线性混合模型")
#bp_data <- read.csv("ChrVNE.2Groups.csv",header = T)
bp_data <- read.csv("ChrVNE.csv",header = T)

mod0 <- lme(VNE ~ group, random = ~1 | id, data = bp_data) #AIC 339.6334

mod1 <- lme(VNE ~ group2, random = ~1 | length, data = bp_data) #AIC 339.6334

mod2 <- lme(VNE ~ group, random = ~1 | species, data = bp_data) #AIC 242.9791 

mod3 <- lme(VNE ~ group + length, random = ~1 | species, data = bp_data) #AIC 150.9893

mod4 <- update(mod3, .~. + group:length) #AIC 228.1895

mod5 <- lme(VNE ~ length, random = ~1 | species, data = bp_data) #AIC 144.5144  但是和mod3没有显著差异

mod6 <- lme(VNE ~ species, random = ~1 | length, data = bp_data)

mod7 <- lme(VNE ~ group + length + species, random = ~1 | id, data = bp_data) #AIC 147.1178

mod8 <- lme(VNE ~ group + species, random = ~1 | id, data = bp_data) #AIC 239.1492



###
bp_data <- read.csv("ChrVNE.csv",header = T)
mod0 <- lmer(VNE~group+(1|species),data = bp_data) #AIC 240

mod1 <- lmer(VNE~group+length+(1|species),data = bp_data) #AIC 105.43

mod2 <- lmer(VNE~length+(1|species),data = bp_data) #AIC 104.25

mod3 <- lmer(VNE~length + species +(1|group),data = bp_data) #AIC 85.301
#最优模型
