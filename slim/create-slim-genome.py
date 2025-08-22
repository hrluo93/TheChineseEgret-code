#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
生成 SLiM initializeGenomicElement(...) 片段（不学习模式；1-based 闭区间）

- g3（间区）可切换两种采样模式：
  1) uniform：长度 ~ Uniform[g3min, g3max]
  2) lognorm：长度 ~ LogNormal(μ, σ)；其中 μ = ln(g3_mean) - 0.5·g3_sigma^2

- g1（外显子）：长度 ~ LogNormal(μ_exon, g1_sigma)，μ_exon = ln(g1_mean) - 0.5·g1_sigma^2
- g2（内含子）：长度 ~ LogNormal(μ_intron, g2_sigma)，μ_intron = ln(g2_mean) - 0.5·g2_sigma^2
- 生成流程：g3 -> g1 -> 循环若干次 (g2, g1)（以概率 0.8 继续）-> 最后补一个 g3
- 所有区段均按 1-based 闭区间输出，即长度 = end - start + 1

注意：
- 本脚本把 -l/--length 解释为**终止坐标**（不是长度）。若从 1 开始且 -l=1_100_000_000，
  则坐标范围为 [1, 1100000000]。
"""

import argparse
import random
import sys
import numpy as np

def simulate_regions(start, L,
                     g3_type, g3min, g3max, g3_mean, g3_sigma,
                     g1_mean, g1_sigma,
                     g2_mean, g2_sigma,
                     p_continue=0.8,
                     min_g1=1, min_g2=10, min_g3=1):
    """
    参数
    ----
    start:    起始坐标（1-based）
    L:        终止坐标（1-based，包含）
    g3_type:  'uniform' 或 'lognorm'
    g3min,g3max: uniform 模式下 g3 的最小/最大长度
    g3_mean,g3_sigma: lognorm 模式下 g3 的原始尺度均值与 log 空间 σ（仅在 g3_type=lognorm 时使用）
    g1_mean,g1_sigma: 外显子 g1 的原始均值与 log 空间 σ
    g2_mean,g2_sigma: 内含子 g2 的原始均值与 log 空间 σ
    p_continue:       是否继续生成 (g2,g1) 对的概率（默认 0.8）
    min_g1,min_g2,min_g3: 采样后对长度的下限修正（避免出现 0 或过短）

    返回
    ----
    [(feature, start, end), ...]
    """
    # 参数检查
    if g1_mean <= 0 or g2_mean <= 0:
        raise ValueError("g1_mean / g2_mean 必须为正数（原始尺度均值）。")
    if g1_sigma <= 0 or g2_sigma <= 0:
        raise ValueError("g1_sigma / g2_sigma 必须为正数（log 空间 σ）。")
    if g3_type == "lognorm":
        if (g3_mean is None) or (g3_sigma is None):
            raise ValueError("g3_type=lognorm 时需要提供 --g3-mean 与 --g3-sigma。")
        if g3_mean <= 0 or g3_sigma <= 0:
            raise ValueError("g3_mean 必须>0 且 g3_sigma 必须>0。")
        mu_g3 = np.log(g3_mean) - 0.5 * (g3_sigma ** 2)

    # 将原始尺度均值换算为对数空间均值 μ： μ = ln(m) - 0.5·σ^2
    mu_exon   = np.log(g1_mean) - 0.5 * (g1_sigma ** 2)
    mu_intron = np.log(g2_mean) - 0.5 * (g2_sigma ** 2)

    # g3 采样器
    if g3_type == "uniform":
        sample_g3 = lambda: random.randint(g3min, g3max)
    else:  # lognorm
        sample_g3 = lambda: max(int(np.random.lognormal(mu_g3, g3_sigma)), min_g3)

    # g1 / g2 采样器
    sample_g1 = lambda: max(int(np.random.lognormal(mu_exon,  g1_sigma)), min_g1)
    sample_g2 = lambda: max(int(np.random.lognormal(mu_intron, g2_sigma)), min_g2)

    segments = []
    base = start

    while base < L:
        # 非编码区 g3
        g3_len = sample_g3()
        if base + g3_len - 1 > L:
            g3_len = L - base + 1
        segments.append(("g3", base, base + g3_len - 1))
        base += g3_len
        if base >= L:
            break

        # 第一个外显子 g1
        ex_len = sample_g1()
        if base + ex_len - 1 > L:
            break
        segments.append(("g1", base, base + ex_len - 1))
        base += ex_len

        # 循环生成 (g2, g1) 对
        while base < L and (random.random() < p_continue):
            in_len = sample_g2()
            if base + in_len - 1 > L:
                break
            segments.append(("g2", base, base + in_len - 1))
            base += in_len
            if base >= L:
                break

            ex_len = sample_g1()
            if base + ex_len - 1 > L:
                break
            segments.append(("g1", base, base + ex_len - 1))
            base += ex_len

    # 尾部补一段 g3（如果还有空间）
    if base < L:
        g3_len = sample_g3()
        if base + g3_len - 1 > L:
            g3_len = L - base + 1
        segments.append(("g3", base, base + g3_len - 1))

    return segments

def output_simulation(sim_segments, out_file):
    lines = [f"initializeGenomicElement({feat}, {s}, {e});"
             for (feat, s, e) in sim_segments]
    if out_file:
        with open(out_file, "w") as f:
            f.write("\n".join(lines) + "\n")
    else:
        print("\n".join(lines))

def main():
    parser = argparse.ArgumentParser(
        description="模拟生成 SLiM 的 initializeGenomicElement() 区段（不学习模式；1-based 闭区间）。"
    )
    # g3 相关
    parser.add_argument("--g3-type", choices=["uniform","lognorm"], default="uniform",
                        help="间区 g3 的长度生成模式：uniform 或 lognorm（默认 uniform）")
    parser.add_argument("-g3min", type=int, default=100,
                        help="g3 uniform 模式最小长度（默认 100）")
    parser.add_argument("-g3max", type=int, default=5000,
                        help="g3 uniform 模式最大长度（默认 5000）")
    parser.add_argument("--g3-mean", type=float,
                        help="g3 lognorm 模式的原始尺度均值（可选）")
    parser.add_argument("--g3-sigma", type=float,
                        help="g3 lognorm 模式的 log 空间 σ（可选）")

    # g1/g2 相关
    parser.add_argument("--g1-mean", type=float, default=50.0,
                        help="外显子 g1 的原始均值 E[X]（默认 50）")
    parser.add_argument("--g1-sigma", type=float, default=0.69,
                        help="外显子 g1 的 log 空间 σ（默认 0.69）")
    parser.add_argument("--g2-mean", type=float, default=100.0,
                        help="内含子 g2 的原始均值 E[X]（默认 100）")
    parser.add_argument("--g2-sigma", type=float, default=0.80,
                        help="内含子 g2 的 log 空间 σ（默认 0.80）")

    # 生成范围与输出
    parser.add_argument("-l", "--length", type=int, required=True,
                        help="终止坐标（例如 1100000000）")
    parser.add_argument("-stat", type=int, default=1,
                        help="起始坐标（1-based，默认 1）")
    parser.add_argument("-o", "--output", default=None,
                        help="输出文件名；未指定则打印到标准输出")

    # 可复现性
    parser.add_argument("--seed", type=int,
                        help="随机种子（可选）")

    args = parser.parse_args()

    # 随机种子
    if args.seed is not None:
        np.random.seed(args.seed)
        random.seed(args.seed)

    # g3 lognorm 参数检查
    if args.g3_type == "lognorm":
        if (args.g3_mean is None) or (args.g3_sigma is None):
            sys.exit("错误：--g3-type lognorm 需要同时提供 --g3-mean 与 --g3-sigma。")

    segs = simulate_regions(
        args.stat, args.length,
        args.g3_type, args.g3min, args.g3max, args.g3_mean, args.g3_sigma,
        args.g1_mean, args.g1_sigma,
        args.g2_mean, args.g2_sigma
    )
    output_simulation(segs, args.output)

if __name__ == "__main__":
    main()
