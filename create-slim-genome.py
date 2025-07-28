#!/usr/bin/env python3
import argparse
import random
import sys
import numpy as np

def simulate_regions(start, L, g3min, g3max,
                     g1_mean, g1_sigma,
                     g2_mean, g2_sigma):
    """
    模拟生成区段：
      - 非编码区 (g3)：长度 ~ Uniform[g3min, g3max]
      - 外显子 (g1)：长度 ~ Lognormal(mu_exon, g1_sigma)，其中
            mu_exon = ln(g1_mean) - 0.5 * g1_sigma^2
        采样后取整并 +1，避免 0 长度。
      - 内含子 (g2)：长度 ~ Lognormal(mu_intron, g2_sigma)，其中
            mu_intron = ln(g2_mean) - 0.5 * g2_sigma^2
        采样后取整并 +10，避免太短。
      - 生成流程：g3 -> g1 -> 循环若干次 (g2, g1)，循环概率 0.8 -> 末尾再加一段 g3。

    参数：
      start: 起始坐标（1-based）
      L: 染色体模拟总长度（终止坐标）
      g3min, g3max: 非编码区 g3 的最小/最大长度（整数）
      g1_mean: 外显子 g1 的目标平均长度（原始尺度 E[X]）
      g1_sigma: 外显子 g1 的对数标准差（log-space σ）
      g2_mean: 内含子 g2 的目标平均长度（原始尺度 E[X]）
      g2_sigma: 内含子 g2 的对数标准差（log-space σ）

    返回：
      列表 [(feature, start, end), ...]，用于 SLiM 的 initializeGenomicElement().
    """
    if g1_mean <= 0 or g2_mean <= 0:
        raise ValueError("g1_mean / g2_mean 必须为正数（原始尺度平均长度）。")
    if g1_sigma <= 0 or g2_sigma <= 0:
        raise ValueError("g1_sigma / g2_sigma 必须为正数（log-space 标准差）。")

    # 将原始尺度均值 m 换算为对数空间均值 μ： μ = ln(m) - 0.5*σ^2
    mu_exon   = np.log(g1_mean) - 0.5 * (g1_sigma ** 2)
    mu_intron = np.log(g2_mean) - 0.5 * (g2_sigma ** 2)

    segments = []
    base = start
    while base < L:
        # 非编码区 g3
        nc_length = random.randint(g3min, g3max)
        if base + nc_length - 1 > L:
            nc_length = L - base + 1
        segments.append(("g3", base, base + nc_length - 1))
        base += nc_length
        if base >= L:
            break

        # 第一个外显子 g1
        ex_length = int(np.random.lognormal(mu_exon, g1_sigma)) + 1
        if base + ex_length - 1 > L:
            break
        segments.append(("g1", base, base + ex_length - 1))
        base += ex_length

        # 循环生成 (g2, g1)
        while base < L and random.random() < 0.8:
            in_length = int(np.random.lognormal(mu_intron, g2_sigma)) + 10
            if base + in_length - 1 > L:
                break
            segments.append(("g2", base, base + in_length - 1))
            base += in_length

            ex_length = int(np.random.lognormal(mu_exon, g1_sigma)) + 1
            if base + ex_length - 1 > L:
                break
            segments.append(("g1", base, base + ex_length - 1))
            base += ex_length

    # 尾部补一段 g3
    if base < L:
        nc_length = random.randint(g3min, g3max)
        if base + nc_length - 1 > L:
            nc_length = L - base + 1
        segments.append(("g3", base, base + nc_length - 1))

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
        description="模拟生成 SLiM 的 initializeGenomicElement() 区段（不学习模式）。"
    )
    parser.add_argument("-g3min", type=int, default=100,
                        help="非编码区(g3)最小长度，默认 100")
    parser.add_argument("-g3max", type=int, default=5000,
                        help="非编码区(g3)最大长度，默认 5000")
    parser.add_argument("-l", "--length", type=int, required=True,
                        help="染色体模拟总长度（终止坐标，例如 1100000000）")
    parser.add_argument("-stat", type=int, default=1,
                        help="起始坐标（1-based），默认 1")
    parser.add_argument("--g1-mean", type=float, default=50.0,
                        help="外显子 g1 的目标平均长度（原始尺度 E[X]），默认 50")
    parser.add_argument("--g1-sigma", type=float, default=0.69,
                        help="外显子 g1 的对数标准差（log-space σ），默认 0.69")
    parser.add_argument("--g2-mean", type=float, default=100.0,
                        help="内含子 g2 的目标平均长度（原始尺度 E[X]），默认 100")
    parser.add_argument("--g2-sigma", type=float, default=0.80,
                        help="内含子 g2 的对数标准差（log-space σ），默认 0.80")
    parser.add_argument("-o", "--output", default=None,
                        help="输出文件名；未指定则打印到标准输出")
    args = parser.parse_args()

    segs = simulate_regions(
        args.stat, args.length,
        args.g3min, args.g3max,
        args.g1_mean, args.g1_sigma,
        args.g2_mean, args.g2_sigma
    )
    output_simulation(segs, args.output)

if __name__ == "__main__":
    main()
if __name__ == "__main__":
    main()
