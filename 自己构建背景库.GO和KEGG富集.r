library(clusterProfiler)
library(dplyr)
library(stringr)
library(enrichplot)#GO,KEGG,GSEA
options(stringsAsFactors = F)
##===============================STEP1:GO注释生成=======================
#自己构建的话，首先需要读入文件
setwd("~")
setwd("E:/workspace/spalax/6.CarmeliFussion/compartment")
setwd("E:/workspace/spalax/4.GaliliFussion/3D.Genomics/compartment/231126_fusion")
egg <- read.delim("galili.eggnog.annotation",header = T,sep="\t")
egg[egg==""]<-NA  #将空行变成NA，方便下面的去除
#从文件中挑出基因query与eggnog注释信息
##gene_info <- egg %>% 
##  dplyr::select(GID = query, GENENAME = eggNOG_OGs) %>% na.omit()
#挑出query_name与GO注释信息
gterms <- egg %>%
  dplyr::select(query, GOs) %>% na.omit()
gene_ids <- egg$query
eggnog_lines_with_go <- egg$GOs!= ""
eggnog_lines_with_go
eggnog_annoations_go <- str_split(egg[eggnog_lines_with_go,]$GOs, ",")
gene2go <- data.frame(gene = rep(gene_ids[eggnog_lines_with_go],
                                 times = sapply(eggnog_annoations_go, length)),
                      term = unlist(eggnog_annoations_go))
names(gene2go) <- c('gene_id', 'ID')
go2name <- read.delim('GO.library', header = FALSE, stringsAsFactors = FALSE)
names(go2name) <- c('ID', 'Description', 'Ontology')

go_anno <- merge(gene2go, go2name, by = 'ID', all.x = TRUE)

## 将GO注释信息保存
#save(go_anno,file = "spalax2Rabbit.NoJac.InputGene_GO.rda")

##===============================STEP3:GO富集分析=================
#目标基因列表(全部基因)
gene_select <- read.delim(file = 'BtoA.region.gene', stringsAsFactors = FALSE,header = F)$V1
#GO 富集分析
#默认以所有注释到 GO 的基因为背景集，也可通过 universe 参数输入背景集
#默认以 p<0.05 为标准，Benjamini 方法校正 p 值，q 值阈值 0.2
#默认输出 top500 富集结果
#如果想输出所有富集结果（不考虑 p 值阈值等），将 p、q 等值设置为 1 即可
#或者直接在 enrichResult 类对象中直接提取需要的结果
go_rich <- enricher(gene = gene_select,
                    TERM2GENE = go_anno[c('ID', 'gene_id')], 
                    TERM2NAME = go_anno[c('ID', 'Description')], 
                    pvalueCutoff = 0.05, 
                    #pAdjustMethod = 'BH', 
                    pAdjustMethod = 'none',
                    qvalueCutoff = 0.5
)
#输出默认结果，即根据上述 p 值等阈值筛选后的
write.csv(go_rich,"galili.Chr1.BtoA.GO.自己建库.csv")
dotplot(go_rich,color = "pvalue")
dotplot(go_rich)
tmp <- merge(go_rich, go2name[c('ID', 'Ontology')], by = 'ID')
tmp <- tmp[c(10, 1:9)]
tmp <- tmp[order(tmp$pvalue), ]

dotplot(tmp,color = "pvalue")
write.table(go_rich, 'tmp.csv', sep = ',', row.names = FALSE, quote = FALSE)

plotGOgraph(go_rich)
#ontology should be one of 'BP', 'MF' or 'CC'

GO = go_rich
#simplyfy函数进行瘦身
#只接受enrichGO
#GO_simple = simplify(GO, cutoff=0.5,by="p.adjust",select_fun=min)
#barplot(GO_simple)
#write.csv(GO_simple,"galili_expand.simplify.csv")

#树状图瘦身
GO2 = pairwise_termsim(GO)
treeplot(GO2)

#网状图瘦身
GO2 <- pairwise_termsim(GO)
enrichplot::emapplot(GO2,showCategory = 25, color = "pvalue", layout = "kk",)
enrichplot::emapplot(GO2,showCategory = 25, color = "pvalue", layout = "kk",clusterFunction = stats::kmeans,group_legend =T,group_category=T)
enrichplot::emapplot(GO2,showCategory = 25, color = "pvalue", layout = "kk",
                     group_legend =F,group_category=F)

#showCategory ：设置展示几个通路及相关差异基因，也可以为字符串，指定感兴趣的通路
#node_label：设置是否展示通路名/基因名 （“category,” “gene,” “all” and “none”)
#layout：设置排布样式，e.g. 'star', 'circle', 'gem', 'dh', 'graphopt', 'grid', 'mds', 'randomly', 'fr', 'kk', 'drl' or 'lgl'

#词云瘦身
#install.packages("wordcloud")
library(wordcloud)
wcdf <- read.table(text = GO$GeneRatio, sep = "/")[1]
wcdf$term <-  GO[,2]
wordcloud(words = wcdf$term, freq = wcdf$V1, scale=(c(4, .1)), colors=brewer.pal(8, "Dark2"), max.words = 25)

#simplyfyEnrichment
library(simplifyEnrichment)
GO_df = as.data.frame(GO)
go_id = GO_df$ID
mat = GO_similarity(go_id, ont = 'CC')

pdf("BtoA.fusion.CC.simclust.pdf",width = 7,height = 9)
df = simplifyGO(mat)
dev.off()


##===============================STEP1:KEGG library=======================
###构建KEGG数据库 - KEGG library
# 需要下载 json文件(这是是经常更新的)
# https://www.genome.jp/kegg-bin/get_htext?ko00001
kegg <- function(json = "ko00001.json") {
  pathway2name <- tibble(Pathway = character(), Name = character())
  ko2pathway <- tibble(Ko = character(), Pathway = character())
  
  kegg <- fromJSON(json)
  
  for (a in seq_along(kegg[["children"]][["children"]])) {
    A <- kegg[["children"]][["name"]][[a]]
    
    for (b in seq_along(kegg[["children"]][["children"]][[a]][["children"]])) {
      B <- kegg[["children"]][["children"]][[a]][["name"]][[b]] 
      
      for (c in seq_along(kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]])) {
        pathway_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["name"]][[c]]
        
        pathway_id <- str_match(pathway_info, "ko[0-9]{5}")[1]
        pathway_name <- str_replace(pathway_info, " \\[PATH:ko[0-9]{5}\\]", "") %>% str_replace("[0-9]{5} ", "")
        pathway2name <- rbind(pathway2name, tibble(Pathway = pathway_id, Name = pathway_name))
        
        kos_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]][[c]][["name"]]
        
        kos <- str_match(kos_info, "K[0-9]*")[,1]
        
        ko2pathway <- rbind(ko2pathway, tibble(Ko = kos, Pathway = rep(pathway_id, length(kos))))
      }
    }
  }
  colnames(ko2pathway) <- c("KO","Pathway")
  save(pathway2name, ko2pathway, file = "kegg_info.RData")
  write.table(pathway2name,"KEGG.library",sep="\t",row.names = F)
}


kegg(json = "ko00001.json")


###上面的代码有漏洞，所以运行一下下面这段，但是本质上还是在构建KEGG library
pathway2name <- tibble(Pathway = character(), Name = character())
ko2pathway <- tibble(Ko = character(), Pathway = character())
KEGG_info='~/ko00001.json'
kegg <- fromJSON(KEGG_info)
for (a in seq_along(kegg[["children"]][["children"]])) {
  A <- kegg[["children"]][["name"]][[a]]
  for (b in seq_along(kegg[["children"]][["children"]][[a]][["children"]])) {
    B <- kegg[["children"]][["children"]][[a]][["name"]][[b]] 
    for (c in seq_along(kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]])) {
      pathway_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["name"]][[c]]
      pathway_id <- str_match(pathway_info, "ko[0-9]{5}")[1]
      pathway_name <- str_replace(pathway_info, " \\[PATH:ko[0-9]{5}\\]", "") %>% str_replace("[0-9]{5} ", "")
      pathway2name <- rbind(pathway2name, tibble(Pathway = pathway_id, Name = pathway_name))
      kos_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]][[c]][["name"]]
      kos <- str_match(kos_info, "K[0-9]*")[,1]
      ko2pathway <- rbind(ko2pathway, tibble(KO = kos, Pathway = rep(pathway_id, length(kos))))
    }
  }
}
###
##===============================STEP2:KEGG注释生成=======================
egg <- read.delim("golani.eggnog.annotation",header = T,sep="\t")
egg[egg==""]<-NA#将空行变成NA以便去除

gene2ko <- egg %>%
  dplyr::select(GID = query, KO = KEGG_ko) %>%
  na.omit()
pathway2name <- read.delim("KEGG.library")
colnames(pathway2name)<-c("Pathway","Name")
gene2ko$KO <- str_replace(gene2ko$KO, "ko:","")
gene2pathway <- gene2ko %>% left_join(ko2pathway, by = "KO",) %>% 
  dplyr::select(GID, Pathway) %>%
  na.omit()
kegg_anno<- merge(gene2pathway,pathway2name,by = 'Pathway', all.x = TRUE)[,c(2,1,3)]
colnames(kegg_anno) <- c('gene_id','pathway_id','pathway_description')
save(kegg_anno,file = "judaei_KEGG.rda")


##===============================STEP4:KEGG注释=================
gene_select <- read.delim('judaei.significant.expand.gene', stringsAsFactors = FALSE,header = F)$V1
#KEGG 富集分析
#默认以所有注释到 KEGG 的基因为背景集，也可通过 universe 参数指定其中的一个子集作为背景集
#默认以 p<0.05 为标准，Benjamini 方法校正 p 值，q 值阈值 0.2
#默认输出 top500 富集结果
kegg_rich <- enricher(gene = gene_select,
                      TERM2GENE = kegg_anno[c('pathway_id', 'gene_id')], 
                      TERM2NAME = kegg_anno[c('pathway_id', 'pathway_description')], 
                      pvalueCutoff = 0.05, 
                      pAdjustMethod = 'BH', 
                      qvalueCutoff = 0.2) 
                      #maxGSSize = 500)
dotplot(kegg_rich)

#输出默认结果，即根据上述 p 值等阈值筛选后的
write.table(kegg_rich, 'judaei.expand.KEGG.csv', sep = ',', row.names = FALSE, quote = FALSE)



#比较富集
setwd("E:/workspace/spalax/4.GaliliFussion/3D.Genomics/compartment/231126_fusion")
b2a <- read.delim(file = 'BtoA.region.gene', stringsAsFactors = FALSE,header = F)$V1
a2b <- read.delim(file = 'AtoB.region.gene', stringsAsFactors = FALSE,header = F)$V1

practice1 = list(b2a, a2b)

xx <- compareCluster(practice1, fun="enricher", 
                     TERM2GENE = go_anno[c('ID', 'gene_id')], 
                     TERM2NAME = go_anno[c('ID', 'Description')], 
                     pvalueCutoff = 0.2, 
                     pAdjustMethod = 'BH', 
                     #pAdjustMethod = 'none',
                     qvalueCutoff = 0.2
                     ) 

dotplot(xx,showCategory=20,includeAll=TRUE) + 
  scale_color_continuous(low='purple', high='green')