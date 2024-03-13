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
install_github("nclark-lab/RERconverge")
help()
install.packages("devtools")
devtools::install_github('lchiffon/REmap')
data("toyTrees.RData")

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
setwd("E:/workspace/spalax/ComparativeGenomics/relax/eye.gene/aln")
in1 = list.files("E:/workspace/spalax/ComparativeGenomics/relax/eye.gene/aln")

estimatePhangornTreeAll(
  alnfiles = in1,
  #alndir = 'E:/workspace/spalax/ComparativeGenomics/relax/eye.gene/aln' ,
  pattern = ".fa.filter",
  "E:/workspace/spalax/ComparativeGenomics/relax/eye.gene/species.tree",
  output.file = "test.out",
  submodel = "LG",
  type = "DNA",
  format = "fasta",
  #k = 4,
)
