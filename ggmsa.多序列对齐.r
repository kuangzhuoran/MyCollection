#devtools::install_github("Ikarto/R4RNA")
#devtools::install_github("YuLab-SMU/ggmsa")

library(ggmsa)
setwd("~")

setwd("E:/workspace/spalax/ComparativeGenomics/relax/BranchLength")

fai<-"all.CL3.fa.mafft.trimal"#读入文

fai<-"test3"#读入文件
ggmsa(fai, seq_name = TRUE)

ggmsa(fai, seq_name = TRUE, show.legend = F,
      consensus_views = T)

ggmsa(fai, border = NA, font = NULL,seq_name = T, color = "LETTER")  + facet_msa(field = 100)

ggmsa(fai, border = NA, seq_name = T, color = "Clustal")  + facet_msa(field = 200)

?ggmsa
