#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
calc_garud_h_from_vcf.py

功能说明
--------
本脚本从 phased 的 VCF 或 VCF.GZ 文件中读取二倍体基因型数据，
将每个样本拆分为两条 haplotype，并在滑动窗口内统计 haplotype 频率，
进而计算 Garud's H statistics，包括：

1. H1
2. H2
3. H12
4. H123
5. H2/H1

脚本支持两种窗口定义方式：
1. 固定物理长度窗口，也就是按 bp 滑窗
2. 固定 SNP 数窗口，也就是按位点数滑窗

同时支持一个可选参数，用于决定是否只保留 focal population 中的
segregating sites。默认不强制过滤，即保留窗口中的所有可用位点。
若设置为只保留 segregating sites，则会去掉在当前群体内已经单态的位点。

适用数据
--------
1. 输入必须是 phased VCF
2. 仅支持 diploid genotype
3. 默认只处理二等位 SNP
4. 含缺失或非 phased 的位点会被跳过
5. 建议输入文件为单条染色体，或已按染色体拆分

统计量定义
----------
设一个窗口内不同 haplotype 的频率为
p1 >= p2 >= p3 >= ...

则：

H1   = sum(pi^2)
H2   = H1 - p1^2
H12  = (p1 + p2)^2 + sum_{i>2}(pi^2)
H123 = (p1 + p2 + p3)^2 + sum_{i>3}(pi^2)
H2/H1 = H2 / H1

实现思路
--------
1. 逐位点读取 VCF
2. 将每个位点的 phased GT 转换为所有 haplotype 的等位基因序列
3. 根据窗口模式提取窗口内位点
4. 将窗口内所有位点拼接为完整 haplotype string
5. 统计不同 haplotype string 的频率分布
6. 根据定义直接计算 H1、H2、H12、H123、H2/H1

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
1. 若窗口已固定，则单态位点本身不会改变 H1、H12、H123、H2/H1 的数值，
   但在固定 SNP 数窗口下，会改变窗口中被纳入的实际位点集合
2. 若样本量较小，或窗口中 SNP 数较少，则这些统计量会较离散，解释时应谨慎
3. 若希望更贴近 Garud 等人的原始 H12 扫描思路，通常建议：
   使用 fixed SNP windows，并只保留 focal population 中的 segregating SNP
"""

import argparse
import gzip
from collections import Counter


def open_textfile(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Calculate H1, H2, H12, H123, and H2/H1 from a phased VCF using sliding windows."
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


def parse_gt_to_haps(gt_field):
    gt = gt_field.split(":")[0]

    if "." in gt:
        return None
    if "|" not in gt:
        return None

    alleles = gt.split("|")
    if len(alleles) != 2:
        return None

    return alleles


def is_biallelic_snp(ref, alt):
    if len(ref) != 1:
        return False
    alt_alleles = alt.split(",")
    if len(alt_alleles) != 1:
        return False
    if len(alt_alleles[0]) != 1:
        return False
    return True


def is_segregating(hap_alleles):
    uniq = set(hap_alleles)
    return len(uniq) > 1


def compute_garud_stats(hap_strings):
    n_haps = len(hap_strings)
    if n_haps == 0:
        return (0.0, 0.0, 0.0, 0.0, "NA", 0)

    counts = Counter(hap_strings)
    freqs = sorted((c / n_haps for c in counts.values()), reverse=True)

    h1 = sum(p * p for p in freqs)
    p1 = freqs[0] if len(freqs) >= 1 else 0.0
    p2 = freqs[1] if len(freqs) >= 2 else 0.0
    p3 = freqs[2] if len(freqs) >= 3 else 0.0

    h2 = h1 - p1 * p1

    if len(freqs) >= 2:
        h12 = (p1 + p2) ** 2 + sum(p * p for p in freqs[2:])
    else:
        h12 = h1

    if len(freqs) >= 3:
        h123 = (p1 + p2 + p3) ** 2 + sum(p * p for p in freqs[3:])
    elif len(freqs) == 2:
        h123 = (p1 + p2) ** 2
    else:
        h123 = h1

    h2_h1 = h2 / h1 if h1 > 0 else "NA"
    k = len(freqs)

    return (h1, h2, h12, h123, h2_h1, k)


def load_vcf(vcf_path, segregating_only=False):
    chroms = []
    positions = []
    site_haps = []
    sample_names = None

    with open_textfile(vcf_path) as f:
        for line in f:
            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                fields = line.rstrip("\n").split("\t")
                sample_names = fields[9:]
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue

            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]

            if not is_biallelic_snp(ref, alt):
                continue

            gts = fields[9:]
            hap_alleles = []

            bad_site = False
            for gt_field in gts:
                parsed = parse_gt_to_haps(gt_field)
                if parsed is None:
                    bad_site = True
                    break
                hap_alleles.extend(parsed)

            if bad_site:
                continue

            if segregating_only and (not is_segregating(hap_alleles)):
                continue

            chroms.append(chrom)
            positions.append(pos)
            site_haps.append(tuple(hap_alleles))

    if sample_names is None:
        raise ValueError("VCF header line not found.")

    if len(positions) == 0:
        raise ValueError("No usable sites found after filtering.")

    return sample_names, chroms, positions, site_haps


def write_output_header(out_handle):
    header = [
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
    ]
    out_handle.write("\t".join(header) + "\n")


def run_bp_windows(chroms, positions, site_haps, n_samples, window_size, step_size, min_sites, out_handle):
    n_haps = n_samples * 2
    n_sites_total = len(positions)

    left = 0
    right = 0
    window_index = 0

    min_pos = positions[0]
    max_pos = positions[-1]

    start = min_pos
    while start <= max_pos:
        end = start + window_size - 1

        while left < n_sites_total and positions[left] < start:
            left += 1
        while right < n_sites_total and positions[right] <= end:
            right += 1

        idx_start = left
        idx_end = right
        n_sites = idx_end - idx_start

        if n_sites >= min_sites:
            hap_strings = [""] * n_haps
            for site_idx in range(idx_start, idx_end):
                alleles = site_haps[site_idx]
                for h in range(n_haps):
                    hap_strings[h] += alleles[h]

            h1, h2, h12, h123, h2_h1, k = compute_garud_stats(hap_strings)

            chrom = chroms[idx_start] if n_sites > 0 else "NA"
            mid = (start + end) // 2

            row = [
                chrom,
                str(window_index),
                str(start),
                str(end),
                str(mid),
                str(n_sites),
                str(n_samples),
                str(n_haps),
                str(k),
                f"{h1:.6f}",
                f"{h2:.6f}",
                f"{h12:.6f}",
                f"{h123:.6f}",
                "NA" if h2_h1 == "NA" else f"{h2_h1:.6f}",
            ]
            out_handle.write("\t".join(row) + "\n")

        window_index += 1
        start += step_size


def run_snp_windows(chroms, positions, site_haps, n_samples, window_size, step_size, min_sites, out_handle):
    n_haps = n_samples * 2
    n_sites_total = len(positions)
    window_index = 0

    start_idx = 0
    while start_idx < n_sites_total:
        end_idx = min(start_idx + window_size, n_sites_total)
        n_sites = end_idx - start_idx

        if n_sites >= min_sites:
            hap_strings = [""] * n_haps
            for site_idx in range(start_idx, end_idx):
                alleles = site_haps[site_idx]
                for h in range(n_haps):
                    hap_strings[h] += alleles[h]

            h1, h2, h12, h123, h2_h1, k = compute_garud_stats(hap_strings)

            chrom = chroms[start_idx]
            window_start = positions[start_idx]
            window_end = positions[end_idx - 1]
            mid = (window_start + window_end) // 2

            row = [
                chrom,
                str(window_index),
                str(window_start),
                str(window_end),
                str(mid),
                str(n_sites),
                str(n_samples),
                str(n_haps),
                str(k),
                f"{h1:.6f}",
                f"{h2:.6f}",
                f"{h12:.6f}",
                f"{h123:.6f}",
                "NA" if h2_h1 == "NA" else f"{h2_h1:.6f}",
            ]
            out_handle.write("\t".join(row) + "\n")

        window_index += 1
        start_idx += step_size


def main():
    args = parse_args()

    segregating_only = (args.segregating_only == "yes")

    sample_names, chroms, positions, site_haps = load_vcf(
        args.vcf,
        segregating_only=segregating_only
    )

    n_samples = len(sample_names)

    with open(args.out, "w") as out_handle:
        write_output_header(out_handle)

        if args.window_type == "bp":
            run_bp_windows(
                chroms=chroms,
                positions=positions,
                site_haps=site_haps,
                n_samples=n_samples,
                window_size=args.window,
                step_size=args.step,
                min_sites=args.min_sites,
                out_handle=out_handle
            )
        else:
            run_snp_windows(
                chroms=chroms,
                positions=positions,
                site_haps=site_haps,
                n_samples=n_samples,
                window_size=args.window,
                step_size=args.step,
                min_sites=args.min_sites,
                out_handle=out_handle
            )


if __name__ == "__main__":
    main()
