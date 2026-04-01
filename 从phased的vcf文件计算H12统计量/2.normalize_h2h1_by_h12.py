#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
本脚本基于 Garud and Rosenberg (2015) 的方法，
对 H2/H1 进行按 H12 条件下的标准化。

背景
----
设：
    Z = H2/H1

H2/H1 的可取范围并不是固定不变的，
而是依赖于对应窗口的 H12 值。
因此，不同位点之间原始的 H2/H1 数值不能直接比较，
尤其是在这些位点的 H12 差异较大时。

为了解决这个问题，文献推导了：

1. 给定 H12 时，H2/H1 的精确理论上界：
       Zmax(H12)

2. 一个便于计算的近似上界：
       Ymax(H12)

并进一步定义了两种标准化统计量：

1. 精确标准化：
       Z' = Z / Zmax(H12)

2. 近似标准化：
       Z'' = Z / Ymax(H12)

功能
----
本脚本读取 Garud H 统计量结果表，
要求输入文件至少包含以下两列：

- H12
- H2_H1

然后在原表基础上新增以下列：

- H2_H1_Zmax_exact
    给定 H12 时 H2/H1 的精确理论上界

- H2_H1_Ymax_approx
    给定 H12 时 H2/H1 的近似理论上界

- H2_H1_norm_exact
    精确标准化后的 H2/H1，即 Z'

- H2_H1_norm_approx
    近似标准化后的 H2/H1，即 Z''

输入
----
输入文件应为制表符分隔文本文件，通常来自以下脚本的输出：

- calc_garud_h_from_vcf.py
- calc_garud_h_scikit_allel.py

输出
----
输出文件为原始结果表加上标准化相关列后的新表，
便于后续进行 peak 比较、作图和结果解释。

注意事项
--------
1. 标准化后的数值理论上应位于 0 到 1 之间
2. 若由于浮点数精度误差导致结果略微小于 0 或大于 1，
   可使用 --clip 参数将结果截断到 [0, 1] 区间
3. H2/H1 的解释应结合 H12 一起进行，
   标准化的意义在 H12 较高的窗口中通常更明显
"""


"""
normalize_h2h1_by_h12.py

Read a Garud H statistics output table and add normalized H2/H1 values
following Garud and Rosenberg (2015).

Input table is expected to contain at least:
- H12
- H2_H1

New output columns:
- H2_H1_Zmax_exact
- H2_H1_Ymax_approx
- H2_H1_norm_exact
- H2_H1_norm_approx

Theory:
Z = H2/H1

Exact normalization:
Z' = Z / Zmax(H12)

Approximate normalization:
Z'' = Z / Ymax(H12)

From Garud and Rosenberg (2015):
Zmax(H12) = (4H12 - 3M^2) / (4H12 - 2M^2)
where
I = ceil((1 + sqrt(8H12 + 1)) / (2H12))
M = [2(I - 1) + 2*sqrt((I^2 - I + 2)H12 - (I + 1))] / (I^2 - I + 2)

Approximate upper bound:
Ymax(H12) = (5 - sqrt(8H12 + 1)) / 4
"""

import argparse
import math
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Normalize H2/H1 by H12 using Garud and Rosenberg (2015)."
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help="Input TSV file produced by calc_garud_h_from_vcf.py or calc_garud_h_scikit_allel.py"
    )
    parser.add_argument(
        "-o", "--output", required=True,
        help="Output TSV file with added normalized columns"
    )
    parser.add_argument(
        "--h12-col", default="H12",
        help="Column name for H12. Default: H12"
    )
    parser.add_argument(
        "--h2h1-col", default="H2_H1",
        help="Column name for H2/H1. Default: H2_H1"
    )
    parser.add_argument(
        "--clip", action="store_true",
        help="Clip normalized values to [0,1] to avoid tiny numerical overflow"
    )
    return parser.parse_args()


def exact_zmax(h12: float) -> float:
    """
    Exact upper bound Zmax(H12) from Garud and Rosenberg (2015), Eq. (12) and Eq. (13).
    Valid for 0 < H12 <= 1.
    """
    if pd.isna(h12):
        return float("nan")
    if h12 <= 0 or h12 > 1:
        return float("nan")

    # Special case for H12 == 1
    if math.isclose(h12, 1.0, rel_tol=1e-12, abs_tol=1e-12):
        return 0.5

    I = math.ceil((1.0 + math.sqrt(8.0 * h12 + 1.0)) / (2.0 * h12))

    denom = I * I - I + 2
    inside = denom * h12 - (I + 1)

    # Numerical tolerance
    if inside < 0 and abs(inside) < 1e-12:
        inside = 0.0
    if inside < 0:
        return float("nan")

    M = (2.0 * (I - 1) + 2.0 * math.sqrt(inside)) / denom

    num = 4.0 * h12 - 3.0 * M * M
    den = 4.0 * h12 - 2.0 * M * M

    if den <= 0:
        return float("nan")

    return num / den


def approx_ymax(h12: float) -> float:
    """
    Approximate upper bound Ymax(H12) from Garud and Rosenberg (2015), Eq. (15).
    Valid for 0 < H12 <= 1.
    """
    if pd.isna(h12):
        return float("nan")
    if h12 <= 0 or h12 > 1:
        return float("nan")

    return (5.0 - math.sqrt(8.0 * h12 + 1.0)) / 4.0


def normalize_value(z: float, upper: float, clip: bool = False) -> float:
    if pd.isna(z) or pd.isna(upper):
        return float("nan")
    if upper <= 0:
        return float("nan")

    out = z / upper
    if clip:
        if out < 0:
            out = 0.0
        elif out > 1:
            out = 1.0
    return out


def main():
    args = parse_args()

    df = pd.read_csv(args.input, sep="\t")

    if args.h12_col not in df.columns:
        raise ValueError(f"Column not found: {args.h12_col}")
    if args.h2h1_col not in df.columns:
        raise ValueError(f"Column not found: {args.h2h1_col}")

    # Convert to numeric safely
    df[args.h12_col] = pd.to_numeric(df[args.h12_col], errors="coerce")
    df[args.h2h1_col] = pd.to_numeric(df[args.h2h1_col], errors="coerce")

    df["H2_H1_Zmax_exact"] = df[args.h12_col].apply(exact_zmax)
    df["H2_H1_Ymax_approx"] = df[args.h12_col].apply(approx_ymax)

    df["H2_H1_norm_exact"] = [
        normalize_value(z, u, clip=args.clip)
        for z, u in zip(df[args.h2h1_col], df["H2_H1_Zmax_exact"])
    ]

    df["H2_H1_norm_approx"] = [
        normalize_value(z, u, clip=args.clip)
        for z, u in zip(df[args.h2h1_col], df["H2_H1_Ymax_approx"])
    ]

    df.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
