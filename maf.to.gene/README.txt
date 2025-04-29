当前目录下的脚本使用GPT修改过
原始的脚本在oringinal下，原始的脚本来自文章"10.1126/science.aav6202"

#以下内容也由GPT生成
# Maf-LST-Gene 处理流程

本项目包含4个脚本，帮助你从大规模MAF多序列比对文件中：
- 提取指定物种的位点信息（LST文件）
- 基于注释GFF定位CDS区域
- 并行提取每个基因的物种序列
- 全流程支持多线程加速，适合大规模数据处理

---

## 脚本列表

### 1. 01.convertMaf2List.pl

功能：
- 从大规模MAF比对中提取指定物种的碑基序列，生成LST表格。
- 以参考物种的柱标坐标为基准，列出各query物种在每个位点的碑基。

参数：
- `-i`   输入maf文件
- `-o`   输出lst文件
- `-ref` 参考物种名
- `-query` 其他物种名（逗号分隔或list文件）
- `-tmp` 参考物种缺失的block输出（默认no_ref.maf）
- `-h/--help` 查看帮助

示例：
```
perl 01.convertMaf2List.pl -i 5zokors.maf -o 5zokors.maf.lst -ref Efontanierii -query Ebaileyi,Esmithi,Erufescens,Erothschildi
```

---

### 2. run_convertMaf_parallel.sh

功能：
- 将超大maf文件按block切分。
- 多线程并行调用 `01.convertMaf2List.pl`。
- 自动合并所有子LST为一个最终大表格。

参数：
- 输入maf文件
- 参考物种
- query物种列表
- 输出文件名
- （可选）每份block数（默认5000）
- （可选）CPU并发数（默认8）

示例：
```
bash run_convertMaf_parallel.sh 5zokors.maf Efontanierii Ebaileyi,Esmithi,Erufescens,Erothschildi 5zokors.merged.lst 5000 16
```

---

### 3. 02.lst2gene.pl

功能：
- 根据LST文件和GFF注释，提取每个基因对应的多物种序列。
- 生成每个基因一个fasta文件，支持正负链互补反向。

参数：
- `-lst` 输入LST文件
- `-gff` 输入GFF注释文件
- `-outdir` 输出目录（默认genes/）
- `-h/--help` 查看帮助

示例：
```
perl 02.lst2gene.pl -lst 5zokors.merged.lst -gff Efontanierii.sorted.rename.gff -outdir genes_output
```

---

### 4. run_lst2gene_parallel.sh

功能：
- 多线程并行调用 `02.lst2gene.pl`的每个基因处理。
- 加速提取基因fasta过程，适合上万基因的大数据量。

参数：
- 输入LST文件
- 输入GFF注释文件
- 输出目录
- （可选）CPU并发数（默认8）

示例：
```
bash run_lst2gene_parallel.sh 5zokors.merged.lst Efontanierii.sorted.rename.gff genes_output 32
```

---

## 推荐整体流程

1. 将大maf文件切块 + 并行提取LST：
    ```
    bash run_convertMaf_parallel.sh 5zokors.maf Efontanierii Ebaileyi,Esmithi,Erufescens,Erothschildi 5zokors.merged.lst 5000 16
    ```

2. 将LST按基因映射提取：
    ```
    bash run_lst2gene_parallel.sh 5zokors.merged.lst Efontanierii.sorted.rename.gff genes_output 32
    ```

---

## 注意事项

- 脚本需要bash环境，perl需安装标准模块（无特殊依赖）。
- 输入maf需为标准maf格式，gff需包含Parent=的CDS注释。
- 建议根据服务器负载合理调整CPU数量。
- 并行过程中临时生成若干小文件，流程结束后会自动清理。
- 如遇内存不足，可调小block_per_file数量。

---

