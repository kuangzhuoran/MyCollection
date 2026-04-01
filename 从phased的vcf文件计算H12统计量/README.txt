===================
如果你用到了这个脚本，请引用：
1. H12 与 H2/H1 的原始文献
Garud NR, Messer PW, Buzbas EO, Petrov DA. 2015.
Recent Selective Sweeps in North American Drosophila melanogaster Show Signatures of Soft Sweeps.
PLoS Genetics 11(2): e1005004.

2. 如果用到了H123，请引用
Harris AM, Garud NR, DeGiorgio M. 2018.
Detection and Classification of Hard and Soft Sweeps from Unphased Genotypes by Multilocus Genotype Identity.
Genetics 210(4): 1429-1452.

Garud's H 统计量计算脚本说明
===========================

本目录包含两个用于从 phased VCF 文件中滑动窗口计算 Garud's H 统计量的 Python 脚本：

1. calc_garud_h_from_vcf.py
2. calc_garud_h_from_vcf_scikit_allel.py

二者都可以计算以下统计量：

- H1
- H2
- H12
- H123
- H2/H1

支持以下两种窗口模式：

- 固定 bp 长度窗口
- 固定 SNP 数窗口

同时支持以下可选过滤：

- 是否只保留 focal population 中的 segregating sites

这两个脚本的主要区别在于：

1. calc_garud_h_from_vcf.py
   为自定义实现版本。
   直接解析 phased VCF 中的 GT 字段，将每个样本拆成两条 haplotype，
   在窗口内统计 haplotype 频率，并按定义直接计算 H1、H2、H12、H123、H2/H1。

2. calc_garud_h_from_vcf_scikit_allel.py
   为调用现成 Python 包 scikit-allel 的版本。
   其中 H1、H12、H123、H2/H1 调用 scikit-allel 的接口计算，
   H2 由 H1 和 H2/H1 反推得到。

----------------------------------------------------------------------
一、适用数据
----------------------------------------------------------------------

这两个脚本都适用于以下数据：

1. 输入文件为 phased VCF 或 VCF.GZ
2. 仅支持 diploid genotype
3. 建议输入文件为单条染色体，或已按染色体拆分
4. 建议位点为二等位 SNP
5. 缺失位点建议预先过滤，虽然脚本也会尽量跳过或过滤异常位点

输入 VCF 中的 GT 字段应类似：

0|0
0|1
1|0
1|1

如果 GT 未定相，例如使用 0/1 而不是 0|1，则不适合本脚本。

----------------------------------------------------------------------
二、统计量定义
----------------------------------------------------------------------

设一个窗口内不同 haplotype 的频率为：

p1 >= p2 >= p3 >= ...

则：

H1   = sum(pi^2)
H2   = H1 - p1^2
H12  = (p1 + p2)^2 + sum_{i>2}(pi^2)
H123 = (p1 + p2 + p3)^2 + sum_{i>3}(pi^2)
H2/H1 = H2 / H1

其中：

- H1 用于衡量窗口内 haplotype 纯合度
- H12 和 H123 对 soft sweep 更敏感
- H2/H1 常用于辅助区分 hard sweep 和 soft sweep

----------------------------------------------------------------------
三、关于 segregating sites
----------------------------------------------------------------------

参数：
--segregating-only {yes,no}

含义：
是否只保留当前 focal population 中真正多态的位点。

默认值：
no

说明：

1. 若窗口已经固定，则单态位点本身不会改变 H1、H12、H123、H2/H1 的数值
2. 但在固定 SNP 数窗口模式下，是否保留单态位点会改变窗口中被纳入的实际位点集合
3. 若希望更贴近 Garud 等人的原始 H12 分析思路，通常建议：
   使用固定 SNP 数窗口，并只保留 focal population 中的 segregating SNP

----------------------------------------------------------------------
四、脚本一：calc_garud_h_from_vcf.py
----------------------------------------------------------------------

说明：
这是自定义实现版本，不依赖专门的 Garud H 统计量包。

优点：
1. 逻辑透明
2. 依赖少
3. 更方便理解和修改

可能缺点：
1. 统计部分由脚本自行实现
2. 若后续扩展更多功能，需要手动维护

基本用法：

1. 固定 bp 窗口

python calc_garud_h_from_vcf.py \
  -i test.vcf \
  -o test.bp.txt \
  -w 20000 \
  -s 5000 \
  --min-sites 1

2. 固定 bp 窗口，只保留 segregating sites
#其实--segregating-only无论怎么设置，对计算出来的H12值是不会影响的
#这个参数只会影响窗口：以固定的SNP数量来定窗口

python calc_garud_h_from_vcf.py \
  -i test.vcf \
  -o test.bp.seg.txt \
  --window-type bp \
  -w 20000 \
  -s 5000 \
  --segregating-only yes \
  --min-sites 1

3. 固定 SNP 数窗口

python calc_garud_h_from_vcf.py \
  -i test.vcf \
  -o test.snp.txt \
  --window-type snp \
  -w 400 \
  -s 50 \
  --segregating-only yes \
  --min-sites 400

----------------------------------------------------------------------
五、脚本二：calc_garud_h_from_vcf_scikit_allel.py
----------------------------------------------------------------------

说明：
这是调用 scikit-allel 现成接口的版本。

依赖包：
- scikit-allel
- numpy
- pandas

安装方式示例：

pip install scikit-allel numpy pandas

或者：

conda install -c conda-forge scikit-allel numpy pandas

优点：
1. 核心 H 统计量调用现成包实现
2. 与 scikit-allel 的数据结构兼容

可能缺点：
1. 依赖额外 Python 包
2. 需要注意不同版本包的接口细节

基本用法：

1. 固定 bp 窗口

python calc_garud_h_from_vcf_scikit_allel.py \
  -i test.vcf \
  -o test.bp.v2.txt \
  -w 20000 \
  -s 5000 \
  --min-sites 1

2. 固定 bp 窗口，只保留 segregating sites

python calc_garud_h_from_vcf_scikit_allel.py \
  -i test.vcf \
  -o test.bp.seg.v2.txt \
  --window-type bp \
  -w 20000 \
  -s 5000 \
  --segregating-only yes \
  --min-sites 1

3. 固定 SNP 数窗口

python calc_garud_h_from_vcf_scikit_allel.py \
  -i test.vcf \
  -o test.snp.v2.txt \
  --window-type snp \
  -w 400 \
  -s 50 \
  --segregating-only yes \
  --min-sites 400

----------------------------------------------------------------------
六、参数说明
----------------------------------------------------------------------

-i, --vcf
输入 phased VCF 或 VCF.GZ 文件

-o, --out
输出结果文件，制表符分隔

--window-type
窗口定义方式，可选：
- bp
- snp

默认：
bp

说明：
1. bp 表示按物理长度定义窗口
2. snp 表示按固定 SNP 数定义窗口

-w, --window
窗口大小

说明：
1. 当 --window-type 为 bp 时，单位为 bp
2. 当 --window-type 为 snp 时，表示每个窗口包含的 SNP 数

-s, --step
窗口步长

说明：
1. 当 --window-type 为 bp 时，单位为 bp
2. 当 --window-type 为 snp 时，表示每次向前滑动的 SNP 数

--segregating-only
是否只保留当前群体中的 segregating sites

可选：
- yes
- no

默认：
no

--min-sites
一个窗口中至少需要多少个位点才输出结果

默认：
1

建议：
1. 若窗口中 SNP 很少，统计量会比较离散
2. 可根据数据密度适当设置为 10、20、50 等

----------------------------------------------------------------------
七、输出文件格式
----------------------------------------------------------------------

输出为 tab 分隔文本文件，字段如下：

chrom
window_index
window_start
window_end
window_mid
n_sites
n_samples
n_haplotypes
K
H1
H2
H12
H123
H2_H1

各字段含义如下：

chrom
染色体名称

window_index
窗口编号，从 0 开始

window_start
窗口起始坐标
对于 bp 模式，为物理窗口起点
对于 snp 模式，为窗口首个 SNP 的位置

window_end
窗口终止坐标
对于 bp 模式，为物理窗口终点
对于 snp 模式，为窗口最后一个 SNP 的位置

window_mid
窗口中点坐标

n_sites
该窗口实际纳入计算的位点数

n_samples
样本数

n_haplotypes
haplotype 数，通常等于 2 × n_samples

K
窗口内不同 haplotype 的个数

H1
Garud's H1

H2
Garud's H2

H12
Garud's H12

H123
Garud's H123

H2_H1
Garud's H2/H1

----------------------------------------------------------------------
八、结果解释注意事项
----------------------------------------------------------------------

1. 样本量较小时，H 统计量会较离散
   例如只有 5 个个体，也就是 10 条 haplotype 时，
   H1、H12、H123、H2/H1 的取值会比较跳

2. 窗口中位点太少时，统计量不稳定
   因此建议设置合适的 --min-sites

3. 固定 bp 窗口与固定 SNP 数窗口的结果不可直接等同
   二者适用于不同分析习惯

4. 若希望更接近原始 H12 文献中的思路，通常推荐：
   fixed SNP windows + segregating-only yes

5. H2/H1 没有普适固定阈值
   通常需要结合 H12、经验分布或模拟结果综合解释

----------------------------------------------------------------------
九、建议使用场景
----------------------------------------------------------------------

如果更看重：
- 代码可读性
- 逻辑透明
- 依赖少
建议优先使用：
calc_garud_h_from_vcf.py

如果更看重：
- 使用现成包
- 核心统计量由已有轮子实现
建议优先使用：
calc_garud_h_from_vcf_scikit_allel.py

----------------------------------------------------------------------
十、备注
----------------------------------------------------------------------

1. 两个脚本的窗口定义、输出格式和参数风格尽量保持一致，便于比较结果
2. 若同一输入下两个脚本结果存在细微差异，通常应重点检查：
   - 窗口边界是否完全一致
   - 是否过滤了相同位点
   - 是否处理了最后一个不完整窗口
   - VCF 中是否包含缺失或异常 GT

3. 建议先在一个较小测试文件上运行，确认输出格式和数值逻辑正确后，再批量跑全基因组数据