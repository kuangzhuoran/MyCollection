#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
calc_garud_h_scikit_allel.py

功能说明
--------
本脚本基于 scikit-allel 提供的 Garud's H statistics 接口，
从 phased 的 VCF 或 VCF.GZ 文件中读取二倍体基因型数据，
并在滑动窗口中计算以下统计量：

1. H1
2. H2
3. H12
4. H123
5. H2/H1

其中：
1. H1、H12、H123、H2/H1 由 scikit-allel 的现成函数计算
2. H2 由定义 H2 = H1 * (H2/H1) 反推得到

脚本支持两种窗口模式：
1. 固定 bp 长度窗口
2. 固定 SNP 数窗口

同时支持一个参数，控制是否只保留 focal population 中的
segregating sites。默认不强制过滤。

依赖包
------
本脚本依赖以下 Python 包：

1. scikit-allel
2. numpy
3. pandas

可使用如下命令安装：

pip install scikit-allel numpy pandas

适用数据
--------
1. 输入必须是 phased VCF
2. 仅支持 diploid genotype
3. 默认假设输入位点已经是二等位位点
4. 含缺失位点会在读入后被去除
5. 建议输入文件为单条染色体，或已按染色体拆分

实现思路
--------
1. 使用 allel.read_vcf() 读入 CHROM、POS、GT、samples
2. 将 genotype array 转为 haplotype array
3. 若指定 segregating-only，则仅保留当前群体中的多态位点
4. 若窗口类型为 snp：
   调用 allel.moving_garud_h() 在固定 SNP 数窗口上计算统计量
5. 若窗口类型为 bp：
   先按物理位置切出窗口，再对每个窗口调用 allel.garud_h()
6. 将结果整理为表格输出

统计量定义
----------
设窗口内不同 haplotype 的频率为
p1 >= p2 >= p3 >= ...

则：

H1   = sum(pi^2)
H2   = H1 - p1^2
H12  = (p1 + p2)^2 + sum_{i>2}(pi^2)
H123 = (p1 + p2 + p3)^2 + sum_{i>3}(pi^2)
H2/H1 = H2 / H1

参数说明
--------
--window-type
    bp 或 snp
    bp 表示窗口大小和步长按物理距离定义
    snp 表示窗口大小和步长按 SNP 数定义

--segregating-only
    yes 或 no
    是否只保留当前群体中的 segregating sites
    默认 no

--min-sites
    一个窗口中至少需要多少个位点才输出结果

输出结果
--------
输出为制表符分隔文件，包含以下字段：

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

其中：
1. K 表示窗口内不同 haplotype 的个数
2. n_sites 表示该窗口纳入计算的位点数
3. 当窗口定义为 bp 时，window_start 和 window_end 为窗口物理边界
4. 当窗口定义为 snp 时，window_start 和 window_end 为窗口首尾 SNP 的物理坐标

注意事项
--------
1. 本脚本的核心 H 统计量计算依赖 scikit-allel
2. H2 不是 scikit-allel 直接返回的值，因此由 H1 和 H2/H1 计算得到
3. 若输入包含多条染色体，建议先按染色体拆分后再运行，尤其是在 bp 窗口模式下
4. 若样本量较小或窗口 SNP 数较少，统计量可能较离散
5. 若希望尽量贴近原始 H12 方法，通常建议采用 fixed SNP windows，
   并只保留当前群体中的 segregating SNP
"""


import argparse
import numpy as np
import pandas as pd
import allel


def parse_args():
    parser = argparse.ArgumentParser(
        description="Calculate H1, H2, H12, H123, and H2/H1 from a phased VCF using scikit-allel."
    )

    parser.add_argument("-i", "--vcf", required=True, help="Input phased VCF or VCF.GZ")
    parser.add_argument("-o", "--out", required=True, help="Output TSV")

    parser.add_argument(
        "--window-type",
        choices=["bp", "snp"],
        default="bp",
        help="Window definition: bp for physical windows, snp for fixed-number-of-SNP windows. Default: bp"
    )
    parser.add_argument(
        "-w", "--window",
        type=int,
        required=True,
        help="Window size. Unit depends on --window-type"
    )
    parser.add_argument(
        "-s", "--step",
        type=int,
        required=True,
        help="Step size. Unit depends on --window-type"
    )

    parser.add_argument(
        "--segregating-only",
        choices=["yes", "no"],
        default="no",
        help="Whether to keep only segregating sites in the focal population. Default: no"
    )

    parser.add_argument(
        "--min-sites",
        type=int,
        default=1,
        help="Minimum number of sites in a window required to report statistics. Default: 1"
    )

    return parser.parse_args()


def load_vcf_as_haplotypes(vcf_path, segregating_only=False):
    callset = allel.read_vcf(
        vcf_path,
        fields=["variants/CHROM", "variants/POS", "calldata/GT", "samples"]
    )

    if callset is None:
        raise ValueError("No variants were read from the VCF.")

    chrom = callset["variants/CHROM"]
    pos = callset["variants/POS"]
    gt = callset["calldata/GT"]
    samples = callset["samples"]

    if gt is None:
        raise ValueError("GT field not found in VCF.")

    # GT shape: (n_variants, n_samples, ploidy)
    g = allel.GenotypeArray(gt)

    # require diploid
    if g.shape[2] != 2:
        raise ValueError("Only diploid phased genotypes are supported.")

    # require no missing
    # g.is_missing() returns shape (n_variants, n_samples)
    miss = g.is_missing().any(axis=1)
    if np.any(miss):
        g = g.compress(~miss, axis=0)
        chrom = chrom[~miss]
        pos = pos[~miss]
    

    # haplotypes shape: (n_variants, n_haplotypes)
    h = g.to_haplotypes()

    # keep only biallelic SNP-like coding in genotype matrix
    # here we assume your VCF already contains only biallelic sites as you said
    # still, after subsetting, segregating_only is optional
    if segregating_only:
        # segregating within focal population: at least two allele states among haplotypes
        is_seg = np.any(h != h[:, [0]], axis=1)
        h = h.compress(is_seg, axis=0)
        chrom = chrom[is_seg]
        pos = pos[is_seg]

    if h.shape[0] == 0:
        raise ValueError("No usable sites remain after filtering.")

    return samples, chrom, pos, h


def compute_h2_from_h1_h2h1(h1, h2_h1):
    h2 = h1 * h2_h1
    return h2


def run_snp_windows(chrom, pos, h, n_samples, window_size, step_size, min_sites):
    # moving_garud_h works on haplotype array with windows along variants
    h1, h12, h123, h2_h1 = allel.moving_garud_h(h, size=window_size, step=step_size)

    results = []
    n_haps = h.shape[1]
    n_variants = h.shape[0]

    starts = np.arange(0, n_variants - window_size + 1, step_size, dtype=int)

    for window_index, start_idx in enumerate(starts):
        end_idx = start_idx + window_size
        n_sites = end_idx - start_idx

        if n_sites < min_sites:
            continue

        window_chrom = chrom[start_idx]
        window_start = int(pos[start_idx])
        window_end = int(pos[end_idx - 1])
        window_mid = (window_start + window_end) // 2

        H1 = float(h1[window_index])
        H12 = float(h12[window_index])
        H123 = float(h123[window_index])
        H2_H1 = float(h2_h1[window_index])
        H2 = float(compute_h2_from_h1_h2h1(H1, H2_H1))

        # K = number of distinct haplotypes in this window
        hw = h[start_idx:end_idx, :]
        # transpose to (n_haps, n_sites)
        hap_rows = hw.T
        K = np.unique(hap_rows, axis=0).shape[0]

        results.append([
            window_chrom,
            window_index,
            window_start,
            window_end,
            window_mid,
            n_sites,
            n_samples,
            n_haps,
            K,
            H1,
            H2,
            H12,
            H123,
            H2_H1
        ])

    return results


def run_bp_windows(chrom, pos, h, n_samples, window_size, step_size, min_sites):
    results = []
    n_haps = h.shape[1]

    min_pos = int(pos[0])
    max_pos = int(pos[-1])

    left = 0
    right = 0
    n_sites_total = len(pos)
    window_index = 0

    start = min_pos
    while start <= max_pos:
        end = start + window_size - 1

        while left < n_sites_total and pos[left] < start:
            left += 1
        while right < n_sites_total and pos[right] <= end:
            right += 1

        idx_start = left
        idx_end = right
        n_sites = idx_end - idx_start

        if n_sites >= min_sites:
            hw = h[idx_start:idx_end, :]
            H1, H12, H123, H2_H1 = allel.garud_h(hw)
            H2 = compute_h2_from_h1_h2h1(H1, H2_H1)

            hap_rows = hw.T
            K = np.unique(hap_rows, axis=0).shape[0]

            window_chrom = chrom[idx_start]
            window_mid = (start + end) // 2

            results.append([
                window_chrom,
                window_index,
                int(start),
                int(end),
                int(window_mid),
                int(n_sites),
                int(n_samples),
                int(n_haps),
                int(K),
                float(H1),
                float(H2),
                float(H12),
                float(H123),
                float(H2_H1)
            ])

        window_index += 1
        start += step_size

    return results


def main():
    args = parse_args()
    segregating_only = (args.segregating_only == "yes")

    samples, chrom, pos, h = load_vcf_as_haplotypes(
        args.vcf,
        segregating_only=segregating_only
    )

    n_samples = len(samples)

    if args.window_type == "snp":
        results = run_snp_windows(
            chrom=chrom,
            pos=pos,
            h=h,
            n_samples=n_samples,
            window_size=args.window,
            step_size=args.step,
            min_sites=args.min_sites
        )
    else:
        results = run_bp_windows(
            chrom=chrom,
            pos=pos,
            h=h,
            n_samples=n_samples,
            window_size=args.window,
            step_size=args.step,
            min_sites=args.min_sites
        )

    df = pd.DataFrame(results, columns=[
        "chrom",
        "window_index",
        "window_start",
        "window_end",
        "window_mid",
        "n_sites",
        "n_samples",
        "n_haplotypes",
        "K",
        "H1",
        "H2",
        "H12",
        "H123",
        "H2_H1"
    ])

    df.to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
