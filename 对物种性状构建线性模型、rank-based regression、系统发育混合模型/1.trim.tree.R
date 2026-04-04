##
library(ape)
setwd("E:/workspace/vision/17.trait.model")

## 输入文件
trait_file <- "test.VA.AL.celldensity.csv"
tree_file <- "fixedname.BigBird.All.NewNames.6714Taxa.replacename.nwk"

## 输出文件
out_tree_file <- "subset_tree.nwk"
out_species_file <- "matched_species_in_tree.txt"
out_missing_file <- "species_not_found_in_tree.txt"

## 读取性状文件
dat <- read.csv(trait_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

## 提取物种名，去掉前后空格
species <- trimws(dat$species)

## 去重
species <- unique(species)

## 读取树
tree <- read.tree(tree_file)

## 统一名字格式
## 性状文件里的空格替换成下划线
species_tree_style <- gsub(" ", "_", species)

## 树上的 tip 名也去掉前后空格
tree$tip.label <- trimws(tree$tip.label)

## 找到匹配的物种
matched_species <- intersect(tree$tip.label, species_tree_style)
missing_species <- setdiff(species_tree_style, tree$tip.label)

## 输出匹配和缺失物种
writeLines(matched_species, out_species_file)
writeLines(missing_species, out_missing_file)

## 从树中提取这些物种
subtree <- keep.tip(tree, matched_species)

## 保存新的子树
write.tree(subtree, file = out_tree_file)

## 打印结果
cat("原树物种数:", length(tree$tip.label), "\n")
cat("性状文件中的物种数:", length(species_tree_style), "\n")
cat("成功匹配的物种数:", length(matched_species), "\n")
cat("未在树中找到的物种数:", length(missing_species), "\n")
cat("子树已保存到:", out_tree_file, "\n")