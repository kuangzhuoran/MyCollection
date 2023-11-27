#ClusterProfile
library(openxlsx)#读取.xlsx文件
library(ggplot2)#柱状图和点状图
library(stringr)#基因ID转换
library(enrichplot)#GO,KEGG,GSEA
library(clusterProfiler)#GO,KEGG,GSEA XDXDXD
library(GOplot)#弦图，弦表图，系统聚类图
library(DOSE)
library(ggnewscale)
library(Rgraphviz)
library(topGO)#绘制通路网络图     

#library(ComplexHeatmap)#绘制图例   包没安上
library(circlize)#绘制富集分析圈图
library(forcats)

library(org.Rn.eg.db)   #Rat
library(org.Mm.eg.db)   #mouse
library(org.Hs.eg.db)   #human
#setwd("E:/workspace/spalax/FuncAnno/golani")
setwd("~")
#setwd("E:/workspace/spalax/4.GaliliFussion/3D.Genomics/compartment/231126_fusion")
#https://mp.weixin.qq.com/s/2_flmfk7LhKeGkzzmrXupQ
#练习用的数据,内容是: Epas1,Egln1,Hsf1,Sftpd  ，即4个基因SYMBOL
#要求格式是#一行#，逗号分隔符，eg :    Epas1,Sftpd,Egln
practice1 = read.csv("BtoA.region.symbol.csv",header = F)
GO = enrichGO(practice1,OrgDb = org.Mm.eg.db,keyType = 'SYMBOL', pvalueCutoff = 0.05
              ,ont = "ALL",pAdjustMethod = 'none',qvalueCutoff = 1)
barplot(GO)
dotplot(GO)
GO = enrichGO(practice1,OrgDb = org.Mm.eg.db,keyType = 'SYMBOL',ont = "CC", pvalueCutoff = 0.05)
write.csv(GO,"BtoA.GO.mouse背景库.csv")

plotGOgraph(GO)

#simplyfy函数进行瘦身
#只接受enrichGO
GO_simple = simplify(GO, cutoff=0.5,by="p.adjust",select_fun=min)
barplot(GO_simple)
write.csv(GO_simple,"galili_expand.simplify.csv")

plotGOgraph(GO_simple)

#树状图瘦身
GO2 = pairwise_termsim(GO)
treeplot(GO2)

#网状图瘦身
GO2 <- pairwise_termsim(GO)
enrichplot::emapplot(GO2,showCategory = 25, color = "pvalue", layout = "kk",)
enrichplot::emapplot(GO2,showCategory = 25, color = "pvalue", layout = "kk",clusterFunction = stats::kmeans,group_legend =T,group_category=T)
#showCategory ：设置展示几个通路及相关差异基因，也可以为字符串，指定感兴趣的通路
#node_label：设置是否展示通路名/基因名 （“category,” “gene,” “all” and “none”)
#layout：设置排布样式，e.g. 'star', 'circle', 'gem', 'dh', 'graphopt', 'grid', 'mds', 'randomly', 'fr', 'kk', 'drl' or 'lgl'

#词云瘦身
install.packages("wordcloud")
library(wordcloud)
wcdf <- read.table(text = GO$GeneRatio, sep = "/")[1]
wcdf$term <-  GO[,2]
wordcloud(words = wcdf$term, freq = wcdf$V1, scale=(c(4, .1)), colors=brewer.pal(8, "Dark2"), max.words = 25)

#simplyfyEnrichment
library(simplifyEnrichment)
GO_df = as.data.frame(GO)
go_id = GO_df$ID
mat = GO_similarity(go_id, ont = 'CC')
df = simplifyGO(mat)


setwd("~/2311.swiss.rescue")
practice0 = read.csv("carmeli_bu.txt",header = F)
practice1 = t(practice0)
practice2 <- bitr(practice1,fromType = 'UNIPROT',toType = 'SYMBOL',OrgDb = org.Rn.eg.db)
write.csv(practice2, file = "carmeli.rescue.rat.tsv")

practice2 <- bitr(practice1,fromType = 'UNIPROT',toType = 'SYMBOL',OrgDb = org.Rn.eg.db)
write.csv(practice2, file = "carmeli.rescue.rat.tsv")

#把基因SYMBOL转换为UNIPROT
#practice2 <- bitr(practice1, fromType='UNIPROT', toType='SYMBOL', OrgDb='org.Rn.eg.db', drop = TRUE) #去除空值
#write.csv(practice2, "rescue.rat.csv")

GO = enrichGO(practice1,OrgDb = org.Mm.eg.db,keyType = 'SYMBOL',ont = "ALL", pvalueCutoff = 0.05)
barplot(GO)
dotplot(GO)
#dotplot(GO,color = "pvalue")

write.csv(GO,"rfmix降低的大INV上的基因的GO富集.csv")

#gene ID转换
practice2 <- bitr(practice1,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Mm.eg.db, drop = T)#去除空值
#KEGG富集要求的是ENTREZID,所以把SYMBOL转换为ENTREZID

#rno mmu hsa
KEGG_database <- 'mmu' #KEGG分析指定物种，物种缩写索引表详见http://www.genome.jp/kegg/catalog/org_list.html
#R.utils::setOption( "clusterProfiler.download.method",'auto' ) rno hsam mmu
R.utils::setOption( "clusterProfiler.download.method",'auto' )
KEGGtry = enrichKEGG(practice2$ENTREZID,organism = 'mmu',pvalueCutoff = 0.5)

barplot(KEGGtry,color = "pvalue") #color = "p.adjust"
dotplot(KEGGtry,color = "pvalue")

write.csv(KEGGtry,"busted鉴定到的正选择基因.KEGG富集.csv")

#查看支持的ID转换类型,如果需要进行ID转换
keytypes(org.Mm.eg.db) #mouse
keytypes(org.Rn.eg.db) #rat
keytypes(org.Hs.eg.db) #homo sapiens


#导入基因id/等等各种id,要求格式是#一行#，逗号分隔符，eg :    Epas1,Sftpd,Egln1
#all = read.csv("09.underground.reg.symbol.fuji.txt", header = F)
#all = read.csv("09.M.reg.symbol.fuji.txt", header = F)
all =  read.csv("M.psg.clusterpro.symbol.txt", header = F)

GO = enrichGO(all,OrgDb = org.Hs.eg.db,keyType = 'SYMBOL',ont = "ALL", 
              pvalueCutoff = 0.4, qvalueCutoff = 0.4)
dotplot(GO) #画点图
barplot(GO)

write.csv(GO,"galili.ITS.loop.GO.csv")

bakegg <- bitr(all,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Mm.eg.db, drop = T)#去除空值
#KEGG富集要求的是ENTREZID,所以把SYMBOL转换为ENTREZID

#KEGG_database <- 'mmu' #KEGG分析指定物种，物种缩写索引表详见http://www.genome.jp/kegg/catalog/org_list.html

KEGGtry = enrichKEGG(bakegg$ENTREZID,organism = 'mmu',pvalueCutoff = 1, qvalueCutoff = 1)
dotplot(KEGGtry)

write.csv(KEGGtry,"KEGG.csv")

enrichplot::cnetplot(GO,circular=F, colorEdge = F, showCategory = 30, layout = "kk")


GO2 <- pairwise_termsim(GO)
enrichplot::emapplot(GO2,showCategory = 25, color = "pvalue", layout = "kk",)
enrichplot::emapplot(GO2,showCategory = 25, color = "pvalue", layout = "kk",clusterFunction = stats::kmeans,group_legend =T,group_category=T)
#showCategory ：设置展示几个通路及相关差异基因，也可以为字符串，指定感兴趣的通路
#node_label：设置是否展示通路名/基因名 （“category,” “gene,” “all” and “none”)
#layout：设置排布样式，e.g. 'star', 'circle', 'gem', 'dh', 'graphopt', 'grid', 'mds', 'randomly', 'fr', 'kk', 'drl' or 'lgl'


write.csv(KEGGtry,"09.underground.reg.KEGG.csv")
#把富集结果保存为csv表格
#2020/06/06 KEGGE富集犯病了
#https://mp.weixin.qq.com/s?__biz=MzAxMDkxODM1Ng==&mid=2247514395&idx=1&sn=97fd9130fdf5fb63a55ce23a336f8313&chksm=9b4bf9a0ac3c70b6de5d56c96cdeb311046ac186eb183bdce7e3ac118c71a942f61b6fcae534&mpshare=1&scene=23&srcid=0606GZMKac345yQvyYBevTz9&sharer_sharetime=1654451188412&sharer_shareid=d348cbfcdb38d11c5fc782cd238cc6a7#rd
#解决方法是运行install.packages('R.utils')
#R.utils::setOption( "clusterProfiler.download.method",'auto' )

b = read.csv("GO_select.csv", header = T)

###比较富集
#https://zhuanlan.zhihu.com/p/573321126
a1 = read.csv("Dscore.MpIncrease.OnMp.Chr1.gene2symbol.csv",header = F)
a <- bitr(a1,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Mm.eg.db, drop = T)#去除空值

b1 = read.csv("Dscore.MpIncrease.OnMp.Chr32.gene2symbol.csv",header = F)
b <- bitr(b1,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Mm.eg.db, drop = T)#去除空值

c1 = read.csv("Dscore.MpDecrease.OnMp.Chr1.gene2symbol.csv",header = F)
c <- bitr(c1,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Mm.eg.db, drop = T)#去除空值

practice1 = list(Mp.Increase.Chr1=a$ENTREZID, Mp.Increase.Chr32=b$ENTREZID, 
                 Mp.Decrease.Chr1=c$ENTREZID)

xx <- compareCluster(practice1, fun="enrichGO", pvalueCutoff=0.05, qvalueCutoff = 0.2, 
                     OrgDb='org.Mm.eg.db') 
xx <- compareCluster(cp, fun="enrichKEGG", organism="mmu", 
                     pvalueCutoff=0.01,qvalueCutoff = 0.05) 


practice1 = list(Mp.Increase.Chr1=a1, Mp.Increase.Chr32=b1, 
                 Mp.Decrease.Chr1=c1)
xx <- compareCluster(practice1, fun="enrichGO", pvalueCutoff=0.05, qvalueCutoff = 0.2, 
                     OrgDb='org.Mm.eg.db',keyType = 'SYMBOL') 

xx <- compareCluster(practice1, fun="enrichGO", pvalueCutoff=0.5, 
                     OrgDb='org.Mm.eg.db',keyType = 'SYMBOL') 
dotplot(xx,showCategory=20,includeAll=TRUE) + 
  scale_color_continuous(low='purple', high='green')


#有些时候需要自己挑选GO term
#把结果导出来, 自己用ggplot2画图
data_GO_sim_fil <- xx@compareClusterResult
df_GO <- write.csv(data_GO_sim_fil,"df_GO.csv")

df_GO <- read.csv("df_GO.csv", header = T)
#手动进行一些处理

library(forcats)
df_GO$Description <- as.factor(df_GO$Description)
df_GO$Description <- fct_inorder(df_GO$Description)

ggplot(df_GO, aes(Description,Cluster)) +
  geom_point(aes(fill=pvalue, size=Count), shape=21)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5),
        axis.text = element_text(color = 'black', size = 12))+
  scale_fill_gradient(low="#56C2BB",high="#EE3029")+
  labs(x=NULL,y=NULL)+
  coord_flip()

################workspace
#eggnog-mapper 2 ClusterProfile
#将注释结果下载到本地，手动删除前三行带井号的行，第四行开头的井号去掉，文件末尾带井号的行去掉
library(stringr)
library(dplyr)

#gene2go
###https://zhuanlan.zhihu.com/p/475588763 改用TBtools也可以，详见此网页；
###如果想用TBtools做富集 https://www.jianshu.com/p/db958ebc73b0

egg<-read.table("FA.emapper.annotations",sep="\t",header=T)  #rggnogmapper的注释结果
gene = read.csv("highfst.csv", header = F)  #筛选到的基因名
gene_ids <- egg$query

#提取GO数据
gene_ids <- egg$query
eggnog_annoations_go <- str_split(egg$GOs, ",")
gene_to_go <- data.frame(gene = rep(gene_ids,
                                    times = sapply(eggnog_annoations_go, length)),
                         term = unlist(eggnog_annoations_go))
gene2go <- filter(gene_to_go, term != "-")

#term2gene <- gene_to_go[,c(2,1)]
term2gene <- gene2go[,c(2,1)]
colnames(term2gene)[1] <- "ID"
df1 <- go2term(term2gene$ID) ##根据ID生成GO注释信息
df2 <- go2ont(term2gene$ID) ##根据ID进行分类：BP, CC, MF
colnames(df1)[1] <- "ID"
colnames(df2)[1] <- "ID"
df <- left_join(term2gene, df1, by = "ID")
df3 <- left_join(df, df2, by = "ID")

#df3_BP <- filter(df3, Ontology == "BP")
#df_three <- split(df3, with(df3, Ontology))
colnames(df3)[4] <- "Class"
gid2gene <- df3[, c("ID", "gene", "Class")]
gid2name <- df3[, c("ID", "Term", "Class")]
ego <- enricher(gene, TERM2GENE = gid2gene, TERM2NAME = gid2name, 
                pvalueCutoff = 0.05, qvalueCutoff = 0.2, minGSSize = 10)
dotplot(ego)
result <- as.data.frame(ego) #转化为数据框方便看富集的结果

gid2gene <- split(gid2gene, with(gid2gene, Class)) #拆分成BP，MF，CC三个数据框
gid2name <- split(gid2name, with(gid2name, Class)) #拆分成BP，MF，CC三个数据框
#gid2gene <- df3_BP[, c("ID", "gene")]

# 2. KEGG pathway
#提取KEGG数据,要注意KO pathway与k number的区别，前者才可转化为代谢功能
gene_ko <- egg %>%
  dplyr::select(GID = query, Ko = KEGG_ko)  #这里的KEGG_ko是K number
gene_ko[,2]<-gsub("ko:","",gene_ko[,2]) #替换去除ko:
gene_ko_list <- str_split(gene_ko$Ko, ",")
gene2ko <- data.frame(gene = rep(gene_ids,
                                 times = sapply(gene_ko_list, length)),
                      term = unlist(gene_ko_list))
gene2ko <- filter(gene2ko, term != "-")

gene_pathway <- egg %>%
  dplyr::select(GID = query, Pathway = KEGG_Pathway)  #这里的才是Pathway
gene_pathway_list <- str_split(gene_pathway$Pathway, ",")
gene2pathway <- data.frame(gene = rep(gene_ids,
                                      times = sapply(gene_pathway_list, length)),
                           term = unlist(gene_pathway_list))
gene2pathway <- filter(gene2pathway, term != "-")

term2gene <- gene2pathway[,c(2,1)]
colnames(term2gene)[1] <- "ID"
pathway2name <- ko2name(term2gene$ID) #根据pathway的id输出注释信息
pathway2name <- na.omit(pathway2name)
pathway2name <- unique.data.frame(pathway2name)
colnames(pathway2name)[1] <- "ID"
ko2gene <- term2gene[grep(pattern = "^ko", term2gene$ID),]
#这里的id有的是map开头的需要去除，因为ko2name只能注释ko开头的

kegg <- enricher(gene, TERM2GENE = ko2gene, TERM2NAME = pathway2name,
                 pvalueCutoff = 0.05, qvalueCutoff = 0.2, minGSSize = 10)
dotplot(kegg)
#browseKEGG(kegg, 'ko05034') #在pathway通路图上标记富集到的基因，会链接到KEGG官网
#

#如果不用R处理得到ko2gene GO2gene gene2GO这些
#可以先用TBtools处理，然后自己再手动修改
#Tbtools处理完的gene2GO是可以直接用的，KEGG则需要去除含有map的行 (eg:  gene1  map00100 这种行就要删除)
#上述处理假设都已经完毕
gene = read.csv("highfst.csv", header = F)
#KEGG富集
#term2pathway内容如下         
#ko001120  gene1
#ko000010  gene2
#map000010  gene3
term2gene = read.table("ko2gene.tsv", header = F, sep = '\t')
colnames(term2gene)[1] <- "ID"
pathway2name <- ko2name(term2gene$ID) #根据pathway的id输出注释信息
pathway2name <- na.omit(pathway2name)
pathway2name <- unique.data.frame(pathway2name)
colnames(pathway2name)[1] <- "ID"
ko2gene <- term2gene[grep(pattern = "^ko", term2gene$ID),]

kegg <- enricher(gene, TERM2GENE = ko2gene, TERM2NAME = pathway2name,
                 pvalueCutoff = 0.5, qvalueCutoff = 0.8,)
dotplot(kegg)
write.csv(kegg,"kegg.csv")
#gene2go内容如下
#GO:000130 gene1
#GO:000120 gene2
#然后这样处理：
term2gene <- gene2go[,c(2,1)]  #clusterProfile格式要求GO号在第一列
colnames(term2gene)[1] <- "ID"
df1 <- go2term(term2gene$ID) ##根据ID生成GO注释信息
df2 <- go2ont(term2gene$ID) ##根据ID进行分类：BP, CC, MF
colnames(df1)[1] <- "ID"
colnames(df2)[1] <- "ID"
df <- left_join(term2gene, df1, by = "ID")
df3 <- left_join(df, df2, by = "ID")

#df3_BP <- filter(df3, Ontology == "BP")
#df_three <- split(df3, with(df3, Ontology))
colnames(df3)[4] <- "Class"
gid2gene <- df3[, c("ID", "gene", "Class")]
gid2name <- df3[, c("ID", "Term", "Class")]

#最后富集用到的就是 上述这两个：TERM2GENE = gid2gene, TERM2NAME = gid2name,
