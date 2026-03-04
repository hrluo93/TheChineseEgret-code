#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
生成 SLiM initializeGenomicElement(...) 片段（1-based 闭区间）

- g3（间区）可切换两种采样模式：
  1) uniform：长度 ~ Uniform[g3min, g3max]
  2) lognorm：长度 ~ LogNormal(μ, σ)；其中 μ = ln(g3_mean) - 0.5·g3_sigma^2

- g1（外显子）：长度 ~ LogNormal(μ_exon, g1_sigma)，μ_exon = ln(g1_mean) - 0.5·g1_sigma^2
- g2（内含子）：长度 ~ LogNormal(μ_intron, g2_sigma)，μ_intron = ln(g2_mean) - 0.5·g2_sigma^2

- 生成流程：g3 -> g1 -> 重复 (g2,g1) 若干次（本版由“概率继续”改为“外显子总数 k 抽样”）-> 最后补一个 g3
- 所有区段均按 1-based 闭区间输出，即长度 = end - start + 1

注意：
- 本脚本 -l/--length 为**终止坐标**（不是长度）。
"""

import argparse
import random
import sys
import numpy as np

def sample_exon_count(mean_val: float, dist: str = "geom", sd: float | None = None) -> int:
    """
    抽样本 gene 的外显子总数 k（>=1）。k=1 表示只有一个 g1；k>1 则有 (k-1) 个 (g2,g1) 对。
    dist:
      - geom   : 支持集 {1,2,...} 的几何分布，E[k]=1/p => p=1/mean
      - poisson: k = 1 + Poisson(λ)，令 E[k]=mean => λ=mean-1
      - fixed  : 固定为 round(mean)
      - nbinom : k = 1 + NB(r,p)；给定均值 mean 和目标 sd（若未给 sd，用温和离散度）
    """
    m = max(float(mean_val), 1.0)
    if dist == "geom":
        p = 1.0 / m
        k = int(np.random.geometric(p))
    elif dist == "poisson":
        lam = max(m - 1.0, 0.0)
        k = 1 + int(np.random.poisson(lam))
    elif dist == "fixed":
        k = max(1, int(round(m)))
    else:
        muY = max(m - 1.0, 0.0)                    # Y ~ NB(r,p)，X = 1 + Y
        if sd is not None and sd > 0:
            vY = sd * sd
        else:
            # 温和离散度：Var ≈ μ + μ^2/10（可按需调整）
            vY = muY + (muY * muY) / 10.0 if muY > 0 else 1.0
        # NB 参数：Var(Y)=μ + μ^2/r => r = μ^2/(Var-μ)，p = r/(r+μ)
        denom = max(vY - muY, 1e-9)
        r = max((muY * muY) / denom, 1e-6)
        p = r / (r + muY) if (r + muY) > 0 else 0.5
        k = 1 + int(np.random.negative_binomial(r, p))
    return max(k, 1)

def simulate_regions(start, L,
                     g3_type, g3min, g3max, g3_mean, g3_sigma,
                     g1_mean, g1_sigma,
                     g2_mean, g2_sigma,
                     exons_per_gene_mean=5.0,
                     exons_per_gene_dist="geom",
                     exons_per_gene_sd=None,
                     min_g1=1, min_g2=10, min_g3=1):
    """
    返回: [(feature, start, end), ...]
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

        # 决定本 gene 的外显子总数 k（>=1），然后追加 (k-1) 个 (g2,g1) 对
        k = sample_exon_count(exons_per_gene_mean, exons_per_gene_dist, exons_per_gene_sd)
        repeats = max(k - 1, 0)

        for _ in range(repeats):
            # g2
            in_len = sample_g2()
            if base + in_len - 1 > L:
                # 空间不够，退出 gene 的构建；while 循环将尝试补尾部 g3
                break
            segments.append(("g2", base, base + in_len - 1))
            base += in_len
            if base >= L:
                break

            # g1
            ex_len = sample_g1()
            if base + ex_len - 1 > L:
                break
            segments.append(("g1", base, base + ex_len - 1))
            base += ex_len

        # 继续 while，进入下一轮 g3/gene 或补尾巴

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

    # 外显子数/基因（新增）
    parser.add_argument("--exons-per-gene-mean", type=float, default=5.0,
                        help="目标每基因外显子数的期望（含第一个 g1），例如 9.6")
    parser.add_argument("--exons-per-gene-dist",
                        choices=["geom","poisson","fixed","nbinom"], default="geom",
                        help="外显子数抽样分布：几何/泊松/固定/负二项（默认几何）")
    parser.add_argument("--exons-per-gene-sd", type=float, default=None,
                        help="仅 nbinom 使用：目标标准差，用于控制离散度")

    # 生成范围与输出
    parser.add_argument("-s", "--start", "-stat", dest="start", type=int, default=1,
                        help="起始坐标（1-based，默认 1；兼容旧参数名 -stat）")
    parser.add_argument("-l", "--length", type=int, required=True,
                        help="终止坐标（例如 1100000000）")
    parser.add_argument("-o", "--output", default=None,
                        help="输出文件名；未指定则打印到标准输出")

    # 可复现性
    parser.add_argument("--seed", type=int, help="随机种子（可选）")

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
        args.start, args.length,
        args.g3_type, args.g3min, args.g3max, args.g3_mean, args.g3_sigma,
        args.g1_mean, args.g1_sigma,
        args.g2_mean, args.g2_sigma,
        exons_per_gene_mean=args.exons_per_gene_mean,
        exons_per_gene_dist=args.exons_per_gene_dist,
        exons_per_gene_sd=args.exons_per_gene_sd
    )
    output_simulation(segs, args.output)

if __name__ == "__main__":
    main()
