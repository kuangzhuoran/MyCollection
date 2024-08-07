library(lme4)#混合效应模型的
library(ggplot2)#绘图的
library(ggeffects) #提取混合效应模型预测结果的
library(extrafont)#添加字体的

###############
#最简单的，一种因素的线性混合模型
#参考https://mp.weixin.qq.com/s/aHwd9SSoiF-acIbUAyYumQ
setwd("D:/KZR/note_2021_5_19/线性混合模型")

mydata = read.csv("示例数据.csv",header = T,row.names = 1)

head(mydata)
#一般情况下，为了使效应值更加显著，我们需要标准化一下我们的解释因子，不标准化也可以的
#我们在这里以不标准化为例，如果想标准化运行下边代码就可以
#mydata$herd_intensity <- scale(mydata$herd_intensity, center = TRUE, scale = TRUE)


#构建混合效应模型的表达式
mixed <- lmer(biomass~herd_intensity + (1|block),data = mydata)
#mixed <- lmer(VNE~group + (1|id),data = mydata)
#函数中指定的模型lmer是线性混合效应模型（LMM）
#固定效应是herd_intensity

#mixed <- lmer(biomass~herd_intensity*block + (1|block),data = mydata)
#此时固定效应是herd_intensity
#随机效应为block,即不同的重复/观测值

#查看模型结果
summary(mixed)


#提取模型的预测结果
pred.mm <- ggpredict(mixed, terms = c("herd_intensity"))
pred.mm
summary(pred.mm)

ggplot(pred.mm) + 
  geom_line(aes(x = x, y = predicted)) +        
  geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), 
              fill = "lightgrey", alpha = 0.5) +  
  geom_point(data = mydata,                     
             aes(x = herd_intensity, y = biomass, colour = Tr,size= herd_intensity)) + 
  scale_colour_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF","#F39B7FFF","#91D1C2FF", "#8491B4FF"))+
  scale_x_continuous(limits = c(0,5))+
  labs(x = "放牧强度", y = "生物量", 
       title = "放牧强度增加导致地上生物量降低") + 
  theme_bw()


#######
#######
#复杂一点的
#参考https://mp.weixin.qq.com/s/nn3sd6IERw1GVQHdg9Blwg
#01.数据集的构建

#假设我们有一个医学研究数据集，涉及一个新的治疗方法对患者血压的影响
#研究设计包括两组：一组接受新治疗（试验组），另一组接受常规治疗（对照组）
#为了评估治疗效果，对每位参与者在治疗前（基线）及治疗后的第1个月、第3个月、第6个月进行血压测量。


#因此，我们的数据集包括以下变量：
#id: 参与者的唯一标识符。
#group: 分组变量，表示参与者属于试验组还是对照组。
#time: 随访时间，有四个水平（基线、1个月、3个月、6个月）。
#BP: 血压测量值，是我们关注的临床结局，属于连续变量。


#针对这个数据集，线性混合效应模型 (LMM) 能够处理这种重复测量数据，因为它可以同时考虑固定效应（如干预组别）和随机效应（如受试者间的随机变异）。
#在这个例子中，固定效应是我们试图估计的总体平均效应（例如，治疗效果和时间的影响），而随机效应允许每个受试者的基线血压和随时间的变化趋势有所不同。


#LMM模型拟合的基本步骤：
#拟合随机截距模型：允许每个参与者的血压基线值不同，但假设*血压随时间的变化趋势*对所有人是相同的。
#增加交互效应：考虑治疗效果可能随时间变化，增加时间和治疗组别的交互效应。
#增加随机斜率：考虑到每个受试者的血压变化趋势可能不同，增加一个随机斜率模型，允许血压随时间的变化斜率在受试者间有所不同。


# 设置随机数种子以保证结果的可复现性
set.seed(123)

# 生成模拟数据
n <- 100 # 参与者总数
times <- c('baseline', '1_month', '3_months', '6_months')
groups <- c('control', 'treatment')
num_times <- length(times)
num_groups <- length(groups)

# 创建参与者ID、分组和时间的数据框
id <- rep(1:n, each = num_times)
group <- rep(rep(groups, each = n / num_groups), each = num_times)
time <- rep(times, times = n)

# 模拟血压测量值
# 基线血压
bp_baseline <- rnorm(n, mean = 120, sd = 10)
# 血压随时间的变化（假设对所有人相同，仅为简化示例）
time_effect <- c(0, -1, -3, -5)

# 治疗效应（对于试验组和对照组不同）
group_effect <- ifelse(group == 'treatment', -2, 0)

# 模拟血压测量值
bp <- numeric(length(id))
for (i in 1:length(id)) {
  time_idx <- match(time[i], times)
  bp[i] <- bp_baseline[id[i]] + time_effect[time_idx] + group_effect[i] + rnorm(1, mean = 0, sd = 5)
}

# 创建数据框
bp_data <- data.frame(id, group, time, BP = bp)

# 检查是否有NA值
anyNA(bp_data)
# 查看前几行数据
head(bp_data)

#id即不同的受试者
#group即是否经过了实验处理
#time是不同的测量时间
#BP是具体的血压，也就是测量值


#02基于nlme包的线性混合模型分析
#2.1 拟合随机截距模型
library(nlme)

mod1 <- lme(BP ~ group + time, random = ~1 | id, data = bp_data)
#此时固定效应是group(即是否治疗)、time(治疗的时间)、以及二者的交互
#随机效应为id,即不同的重复/观测值
summary(mod1)

#解读：模型mod1（随机截距模型）
#固定效应（Fixed effects）: 这部分展示了固定效应变量的估计值、标准误、自由度（DF）、t值和p值。
#在我们的模型中，固定效应包括group和time以及它们对血压（BP）的影响。
  #(Intercept): 表示基线血压的平均估计值。
  #grouptreatment: 治疗组相对于对照组的血压差异，其p-值显示这个差异不显著。
  #time: 不同时间点的血压变化。time6_months和timebaseline显著，说明血压随时间的变化是有统计学意义的。

#随机效应（Random effects）: 描述了随机效应的方差成分，包括截距的标准差和残差的标准差。
#这表示受试者间的基线血压存在差异，以及测量误差或个体响应的变异性。

#模型拟合指标（AIC, BIC, logLik）: 用于比较不同模型的拟合优度。AIC和BIC越低表示模型越优



#2.2增加交互效应：：考虑治疗效果可能随时间变化，增加时间和治疗组别的交互效应

# 增加交互效应
mod2 <- update(mod1, .~. + group:time)
summary(mod2)

# 检验显著性
anova(mod2)
#解读：模型mod2（增加交互效应）
#固定效应部分：中新增了交互项，用于评估治疗效果随时间变化的差异。
#交互项的p-值较高，表明交互效应不显著。group:time

#随机效应部分：与mod1相似，表示模型在考虑交互效应后，随机效应的组成没有显著变化。


#2.3增加随机斜率
#考虑到每个受试者的血压变化趋势可能不同，我们增加一个随机斜率模型，允许血压随时间的变化斜率在受试者间有所不同。

# 增加随机斜率
mod3 <- update(mod1, random = ~1 + time | id)
summary(mod3)

#解读：模型mod3（增加随机斜率）
#随机效应中，增加了时间的斜率：允许不同受试者对时间的响应（血压变化趋势）不同。
#这里的标准差和相关系数描述了随机斜率的变异性和不同时间点间的相关性。

#固定效应的解释与前面模型相同：但是注意到随机斜率模型可能会影响固定效应的解释。


#最后，比较mod1和mod3，看看是否有必要加入随机斜率：
# 比较不同模型
anova(mod1, mod3)


#解读：比较不同模型
#使用anova函数比较mod1和mod3时，我们关注以下几点：
#L.Ratio（似然比）: 表示两个模型拟合优度的差异，一个较大的值通常表示较复杂的模型（例如增加了随机斜率的mod3）有更好的拟合。
#p-value: 比较两个模型是否存在显著差异。在这个例子中，p-值为0.2383，表明增加随机斜率的模型（mod3）与只有随机截距的模型（mod1）相比，并没有显著提高模型的拟合度。



