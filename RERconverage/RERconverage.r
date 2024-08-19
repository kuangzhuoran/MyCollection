library(RERconverge)
library(devtools)
library(readTrees)
##检查安装RERconverage了米有
if (!require("RERconverge", character.only=T, quietly=T)) {
  require(devtools)
  install_github("nclark-lab/RERconverge", ref="master")
  #"ref" can be modified to specify a particular branch
}
library(RERconverge)
rerpath = find.package('RERconverge') #If this errors, there is an issue with installation
print(rerpath)
##安装RERconverge
#install_github("nclark-lab/RERconverge")
#help()
#install.packages("devtools")
#devtools::install_github('lchiffon/REmap')
#data("toyTrees.RData")

##读入文件
#?readTrees
toytreefile = "all.txt"
toyTrees=readTrees(toytreefile, max.read = 13287)


setwd("E:/workspace/spalax/ComparativeGenomics/relax/eye.gene")
#toytreefile2 = "subsetMammalGeneTrees.txt"
#toyTrees2=readTrees(toytreefile2, max.read = 20)

data("logAdultWeightcm")

mamRERw = getAllResiduals(toyTrees2,useSpecies=names(logAdultWeightcm),
                          transform = "sqrt", weighted = T, scale = T)

mamRERw = getAllResiduals(toyTrees,useSpecies = c("carmeli","galili","golani","judaei"),
                          transform = "sqrt", weighted = T, scale = T,plot = T,min.sp=1)

saveRDS(mamRERw, file="spalaxRERw.rds")
newmamRERw = readRDS("spalaxRERw.rds")

noneutherians <- c("carmeli","judaei","galili","golani")
XD = c("Rabbit") #指定外群
par(mfrow=c(1,2))
avgtree=plotTreeHighlightBranches(toyTrees$masterTree, outgroup=XD,
                                  hlspecies=c("Vole","Squirrel"), hlcols=c("blue","red"),
                                  main="Average tree") #plot average tree
bend3tree=plotTreeHighlightBranches(toyTrees$trees$evm.model.ptg000002l.100, outgroup=XD,
                                    hlspecies=c("Vole","Squirrel"), hlcols=c("blue","red"),
                                    main="Trpm1 tree") #plot individual gene tree


par(mfrow=c(1,1))
phenvExample <- foreground2Paths(c("carmeli","galili","golani","judaei"),
                                 toyTrees,clade="ancestral",plotTree = T)

plotRers(mamRERw,"evm.model.ptg000044l.99",phenv=phenvExample,plot=T) #plot RERs

bend3rers = returnRersAsTree(toyTrees, mamRERw, "evm.model.ptg000044l.99", plot = TRUE,
                             phenv=phenvExample) #plot RERs

###
##From Others
setwd("E:/workspace/spalax/ComparativeGenomics/relax/eye.gene/aln")

##RERconverage要求基因树和物种树相同，所以先将基因树拓扑结构强行扭为物种树

estimatePhangornTreeAll(
  alndir = 'E:/workspace/spalax/ComparativeGenomics/relax/eye.gene/aln' ,
  pattern = ".fa.filter",
  treefile = "E:/workspace/spalax/ComparativeGenomics/relax/eye.gene/aln/species.tree",
  output.file = "test.out",
  submodel = "GTR",
  type = "DNA",
)

#estimatePhangornTreeAll(alndir = "test",treefile = "Phylogenetic_tree/Hystricomorpha.trees",output.file = "gene.trees",
#                        type = "DNA",submodel="GTR")

#用readTrees读取基因树
#也就是上一步的输入出文件
toyTrees = readTrees("E:/workspace/spalax/ComparativeGenomics/relax/eye.gene/aln/test.out")

#用getAllResiduals估计相对进化率(RER)
mamRERw = getAllResiduals(toyTrees,
                          #useSpecies = c("carmeli","galili","golani","judaei"),
                          transform = "sqrt", weighted = T, scale = T,plot = T,min.sp=1)
#useSpecies:一个矢量，可用于指定要在分析中使用的物种子集



#如果希望保存这个RER对象以供以后使用，可以使用R的saveRDS函数。这将允许您稍后使用readRDS加载它，如果您愿意，可以使用不同的名称
saveRDS(mamRERw, file="spalaxRERw.rds")
newmamRERw = readRDS("spalaxRERw.rds")


noneutherians <- c("carmeli","judaei","galili","golani") #目标物种
XD = c("Rabbit") #指定外群
par(mfrow=c(1,2))
avgtree=plotTreeHighlightBranches(toyTrees$masterTree, outgroup=XD,
                                  hlspecies=c("Vole","Squirrel"), hlcols=c("blue","red"),
                                  main="Average tree") #plot average tree
bend3tree=plotTreeHighlightBranches(toyTrees$trees$evm.model.ptg000002l.143.fa, outgroup=XD,
                                    hlspecies=c("Vole","Squirrel"), hlcols=c("blue","red"),
                                    main="Trpm1 tree") #plot individual gene tree
#左图是一个树形图，其分支长度代表了所有基因的平均比率。右图是同一棵树，但分支长度代表了特定基因的特异性比率。


par(mfrow=c(1,1))
phenvExample <- foreground2Paths(c("carmeli","judaei","galili","golani"),toyTrees,clade="terminal")

plotRers(mamRERw,"evm.model.ptg000002l.143.fa",phenv=phenvExample,plot=T)
#加了phenv=phenvExample，会多一个rho和p值，也就是foreground2Paths中指定的物种们，进化速率是不是更高

plotRers(mamRERw,"evm.model.ptg000002l.143.fa",plot=T)


#我们还可以使用returnRersAsTree函数保存给定基因的所有RERs。
#如果我们包含plot=TRUE，这个函数还将生成一个基因树的梯形图，其中的分支用它们的RERs标记。

#plot RERs as tree
par(mfrow=c(1,1))
bend3rers = returnRersAsTree(toyTrees, mamRERw, "evm.model.ptg000002l.143.fa", plot = TRUE,
                             phenv=phenvExample) #plot RERs

#此函数还返回一个树，其分支长度表示该基因的RERs。我们可以使用ape包函数write.tree打印或保存树:
strwrap(gsub(":",write.tree(bend3rers),replacement=": "))
write.tree(bend3rers, file='test.rer.out.nwk')


#函数returnRersAsTreesAll将生成一个“multiPhylo”类对象，其中每个元素包含一个命名的基因树，其分支长度表示该基因的RERs。
multirers = returnRersAsTreesAll(toyTrees,mamRERw)
write.tree(multirers, file='toyRERs.nwk', tree.names=TRUE)
#然后，这可以用于进一步分析或自定义绘图，并且可以使用write.tree将其保存到文件中。


#另一个用于可视化给定树的RERs的有用函数是treePlotRers。
#当使用type = label调用时，它将生成与上面所示相同的图，即一个分支由RER标记的梯形图。
#或者，当使用type = color调用时，它将生成一个渐变图，其中的RERs在分支上显示为颜色热图

#visualize RERs along branches as a heatmap
newbend3rers = treePlotRers(treesObj=toyTrees, rermat=mamRERw, index="evm.model.ptg000002l.143.fa",
                            type="c", nlevels=9, figwid=10)
