#!/usr/bin/env python3
import argparse
import random
import sys
import numpy as np

def simulate_regions(start, L, g3min, g3max):
    """
    模拟生成区段：
      - 非编码区 (g3)：长度随机取自均匀分布 [g3min, g3max]；
      - 外显子 (g1)：使用对数正态分布采样，参数：mean=log(50), sigma=log(2)，再取整后加1；
      - 内含子 (g2)：使用对数正态分布采样，参数：mean=log(100), sigma=log(1.5)，取整后加10；
      - 先生成一段 g3，然后生成 g1，再循环生成 (g2, g1) 对，循环概率 0.8，最后生成一段 g3 作为结束区域。
    参数：
      start: 起始坐标（1-based）
      L: 染色体模拟总长度
      g3min, g3max: 非编码区(g3)最小和最大长度
    返回一个列表，元素为 (feature, start, end) 的三元组。
    """
    segments = []
    base = start
    while base < L:
        # 生成非编码区 g3
        nc_length = random.randint(g3min, g3max)
        if base + nc_length - 1 > L:
            nc_length = L - base + 1
        segments.append(("g3", base, base + nc_length - 1))
        base += nc_length
        if base >= L:
            break
        # 生成第一个外显子 (g1)
        ex_length = int(np.random.lognormal(np.log(50), np.log(2))) + 1
        if base + ex_length - 1 > L:
            break
        segments.append(("g1", base, base + ex_length - 1))
        base += ex_length
        # 循环生成内含子–外显子 (g2, g1) 对
        while base < L and random.random() < 0.8:
            in_length = int(np.random.lognormal(np.log(100), np.log(1.5))) + 10
            if base + in_length - 1 > L:
                break
            segments.append(("g2", base, base + in_length - 1))
            base += in_length
            ex_length = int(np.random.lognormal(np.log(50), np.log(2))) + 1
            if base + ex_length - 1 > L:
                break
            segments.append(("g1", base, base + ex_length - 1))
            base += ex_length
    # 最后剩余部分生成一段非编码区
    if base < L:
        nc_length = random.randint(g3min, g3max)
        if base + nc_length - 1 > L:
            nc_length = L - base + 1
        segments.append(("g3", base, base + nc_length - 1))
    return segments

def output_simulation(sim_segments, out_file):
    out_lines = []
    for seg in sim_segments:
        feature, s, e = seg
        out_lines.append(f"initializeGenomicElement({feature}, {s}, {e});")
    if out_file:
        with open(out_file, "w") as f:
            f.write("\n".join(out_lines) + "\n")
    else:
        print("\n".join(out_lines))

def main():
    parser = argparse.ArgumentParser(
        description="模拟生成 SLiM 格式的区段（不学习模式）。\n"
                    "使用 -g3min 和 -g3max 设置非编码区(g3)的最小和最大长度，"
                    "使用 -l 指定染色体总长度，-stat 指定起始坐标，-o 指定输出文件。"
    )
    parser.add_argument("-g3min", type=int, default=100, help="非编码区(g3)最小长度，默认 100")
    parser.add_argument("-g3max", type=int, default=5000, help="非编码区(g3)最大长度，默认 5000")
    parser.add_argument("-l", "--length", type=int, required=True, help="染色体模拟总长度")
    parser.add_argument("-stat", type=int, default=1, help="起始坐标，默认 1")
    parser.add_argument("-o", "--output", default=None, help="输出文件名，若不指定则打印到标准输出")
    
    args = parser.parse_args()
    
    sim_segments = simulate_regions(args.stat, args.length, args.g3min, args.g3max)
    output_simulation(sim_segments, args.output)

if __name__ == "__main__":
    main()
