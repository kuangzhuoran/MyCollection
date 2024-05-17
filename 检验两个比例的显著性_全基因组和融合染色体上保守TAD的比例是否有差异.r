mytable3 <- as.table(rbind(c(11665, 484), c(3478, 199)))
             
dimnames(mytable3) <- list(index = c("0", "1"),
                            Death = c("0","1"))  

chisq.test(mytable3, correct=T)
chisq.test(mytable3, correct=F)


fisher.test(mytable3)


#双比率Z-test
prop.test(n = c(11665, 484), x = c(3478, 199), alternative = 'two.sided', conf.level = 0.95)






#############
#59755是全基因组的保守loop数量，3438是Chr3/融合染色体的保守loop数量
#77506是全基因组的Other loop数量，1175是Chr3/融合染色体的Other loop数量
#prop.test(n = c(59755, 3438), x = c(77506, 1175), alternative = 'two.sided', conf.level = 0.95)

#2781是全基因组保守TAD的数量，97是融合染色体上保守TAD的数量
#3187是全基因组不保守TAD的数量，233是融合染色体上不保守TAD的数量
#mytable3 <- as.table(rbind(c(2781, 97), c(3187, 233)))
mytable3 <- as.table(rbind(c(1346, 38), c(3971, 212)))
                     
dimnames(mytable3) <- list(index = c("0", "1"), Death = c("0","1"))  

chisq.test(mytable3)

fisher.test(mytable3)
                     