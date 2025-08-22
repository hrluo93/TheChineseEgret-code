#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
learn_slim.py

学习模式（-learn）：
  从已有 SLiM initializeGenomicElement(...) 文件中，在指定学习区间内拟合
  g1/g2/g3 的长度与 (g2,g1) 对的数量（基因块大小）与位置的关系（简单 MLP 回归），
  再按 g3→g1→(g2,g1)*→g3 生成新的片段序列，流式写入输出文件。

回退分布（fallback）：
  - g1/g2：对数正态（可通过 --g1-mean/--g1-sigma, --g2-mean/--g2-sigma 提供
            原始尺度均值 + log 空间σ；若未提供，则回到 SLiM 手册示例：
            g1 ~ LogNormal(ln(50), ln(2)), g2 ~ LogNormal(ln(100), ln(1.5))）
  - g3：均匀分布 U[g3min, g3max]

上限保险栅栏（可配置）：
  --max-g1（默认 5000 bp）, --max-g2（默认 200000 bp）, --max-pairs（默认 30）
"""

import argparse
import random
import sys
import re
import numpy as np
import tensorflow as tf
from tensorflow import keras
from keras import layers, models

# -----------------------------
# 解析一行 SLiM 片段
# -----------------------------
def parse_slim_line(line):
    """
    解析类似：
      initializeGenomicElement(g3, 188618852, 188693529);
    返回 (feature, start, end) 或 None
    """
    pattern = r'initializeGenomicElement\(\s*([^,]+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\);'
    m = re.match(pattern, line)
    if m:
        feature = m.group(1).strip()
        start = int(m.group(2))
        end = int(m.group(3))
        return (feature, start, end)
    return None

# -----------------------------
# 读取训练数据（长度）
# -----------------------------
def read_length_data(input_file, inv_interval, stat, total_length):
    """
    对每条记录（g1, g2, g3），若其中点在 inv_interval 内，则记录 (x, length)
    x = (mid - stat)/total_length
    返回 {feature: [(x, length), ...]}
    """
    data = {"g1": [], "g2": [], "g3": []}
    inv_start, inv_end = inv_interval
    with open(input_file) as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            seg = parse_slim_line(line)
            if not seg:
                continue
            feat, s, e = seg
            mid = (s + e) // 2
            if inv_start <= mid <= inv_end and feat in data:
                x = (mid - stat) / float(total_length)
                length = e - s + 1
                data[feat].append((x, length))
    return data

# -----------------------------
# 读取训练数据（基因块对数）
# -----------------------------
def read_count_data(input_file, inv_interval, stat, total_length):
    """
    将 slim 文件按“g3 边界”切成基因块：
      g3 | g1 (g2,g1)* | g3
    在每块内统计 (g2,g1) 的对数（pair_count），用第一段 g1 的中点作为位置代表。
    仅取第一段 g1 中点落在 inv_interval 的块。
    返回 [(x, pair_count), ...]
    """
    blocks = []
    lines = [l.strip() for l in open(input_file) if l.strip()]
    segs = []
    for line in lines:
        p = parse_slim_line(line)
        if p:
            segs.append(p)
    if not segs:
        return blocks

    g3_idx = [i for i, seg in enumerate(segs) if seg[0] == "g3"]
    for i in range(len(g3_idx) - 1):
        block = segs[g3_idx[i]: g3_idx[i+1]]  # 含前不含后
        if len(block) < 3:
            continue
        core = block[1:-1]
        if not core or core[0][0] != "g1":
            continue
        # 统计 (g2, g1) 对
        pair_count = 0
        j = 1
        while j < len(core) - 1:
            if core[j][0] == "g2" and core[j+1][0] == "g1":
                pair_count += 1
                j += 2
            else:
                j += 1
        g1_seg = core[0]
        g1_mid = (g1_seg[1] + g1_seg[2]) // 2
        if inv_interval[0] <= g1_mid <= inv_interval[1]:
            x = (g1_mid - stat) / float(total_length)
            blocks.append((x, pair_count))
    return blocks

# -----------------------------
# 简单回归模型
# -----------------------------
def build_regression_model(input_dim=1):
    model = models.Sequential([
        layers.Input(shape=(input_dim,)),
        layers.Dense(16, activation='relu'),
        layers.Dense(8, activation='relu'),
        layers.Dense(1)
    ])
    model.compile(optimizer='adam', loss='mse')
    return model

def train_length_models(input_file, inv_interval, stat, total_length, epochs=20):
    data = read_length_data(input_file, inv_interval, stat, total_length)
    out = {}
    for feat in ["g1", "g2", "g3"]:
        rows = data.get(feat, [])
        if len(rows) == 0:
            print(f"Warning: No training data for {feat} in learning interval!")
            out[feat] = None
            continue
        arr = np.array(rows, dtype=np.float32)
        X = arr[:, 0:1]
        y = arr[:, 1:2]
        model = build_regression_model()
        model.fit(X, y, epochs=epochs, verbose=0)
        pred = model.predict(X, verbose=0)
        rmse = float(np.sqrt(np.mean((pred.flatten() - y.flatten())**2)))
        print(f"Feature {feat}: trained on {len(X)} samples, RMSE = {rmse:.2f}")
        out[feat] = (model, rmse)
    return out

def train_count_model(input_file, inv_interval, stat, total_length, epochs=20):
    rows = read_count_data(input_file, inv_interval, stat, total_length)
    if len(rows) == 0:
        print("Warning: No count training data in learning interval!")
        return None
    arr = np.array(rows, dtype=np.float32)
    X = arr[:, 0:1]
    y = arr[:, 1:2]
    model = build_regression_model()
    model.fit(X, y, epochs=epochs, verbose=0)
    pred = model.predict(X, verbose=0)
    rmse = float(np.sqrt(np.mean((pred.flatten() - y.flatten())**2)))
    print(f"Count model: trained on {len(X)} samples, RMSE = {rmse:.2f}")
    return (model, rmse)

# -----------------------------
# 预测与回退
# -----------------------------
def predict_length(feature, base, stat, total_length, model_info, fallback_fn, min_val):
    """
    预测长度：若有模型，用模型 + 高斯噪声(rmse)；否则用 fallback_fn()
    """
    if model_info is not None:
        model, rmse = model_info
        x = np.array([[(base - stat) / float(total_length)]], dtype=np.float32)
        pred = float(model.predict(x, verbose=0)[0, 0])
        length = int(round(pred + np.random.normal(0.0, rmse)))
    else:
        length = int(fallback_fn())
    return max(length, min_val)

def predict_count(base, stat, total_length, model_info, fallback_count_fn):
    """
    预测 (g2,g1) 对数：若有模型，用模型 + 噪声；否则回退
    """
    if model_info is not None:
        model, rmse = model_info
        x = np.array([[(base - stat) / float(total_length)]], dtype=np.float32)
        pred = float(model.predict(x, verbose=0)[0, 0])
        cnt = int(round(pred + np.random.normal(0.0, rmse)))
    else:
        cnt = int(fallback_count_fn())
    return max(cnt, 0)

def fallback_count():
    """几何式回退：以 0.8 概率继续叠加对数。"""
    c = 0
    while random.random() < 0.8:
        c += 1
    return c

# -----------------------------
# 流式生成 & 写出
# -----------------------------
def deep_simulate_regions_stream(stat, total_length, length_models, count_model,
                                 g3min, g3max,
                                 g1_mean=None, g1_sigma=None,
                                 g2_mean=None, g2_sigma=None,
                                 fh=None,
                                 max_g1=5000, max_g2=200000, max_pairs=30):
    """
    流式生成并写到 fh：
      while base < total_length:
        g3 -> g1 -> repeat (g2,g1) 'count' times -> 下一块
      末尾若有余量再补一个 g3
    """
    assert fh is not None, "fh (file handle) must be provided."

    # g3 回退
    fb_g3 = lambda: random.randint(int(g3min), int(g3max))

    # g1 回退（对数正态）
    if (g1_mean is not None) and (g1_sigma is not None):
        mu_ex = np.log(float(g1_mean)) - 0.5 * (float(g1_sigma) ** 2)
        fb_g1 = lambda: int(np.random.lognormal(mu_ex, float(g1_sigma))) + 1
    else:
        fb_g1 = lambda: int(np.random.lognormal(np.log(50.0), np.log(2.0))) + 1

    # g2 回退（对数正态）
    if (g2_mean is not None) and (g2_sigma is not None):
        mu_in = np.log(float(g2_mean)) - 0.5 * (float(g2_sigma) ** 2)
        fb_g2 = lambda: int(np.random.lognormal(mu_in, float(g2_sigma))) + 10
    else:
        fb_g2 = lambda: int(np.random.lognormal(np.log(100.0), np.log(1.5))) + 10

    base = int(stat)
    end_limit = int(total_length)

    # 统计信息
    n_g1 = n_g2 = n_g3 = 0
    cut_g1 = cut_g2 = 0
    cut_pairs = 0
    gene_blocks = 0

    def write_seg(tag, s, e):
        fh.write(f"initializeGenomicElement({tag}, {s}, {e});\n")

    while base < end_limit:
        # g3
        g3_len = predict_length("g3", base, stat, total_length,
                                length_models.get("g3"), fb_g3, min_val=g3min)
        if base + g3_len - 1 > end_limit:
            g3_len = end_limit - base + 1
        write_seg("g3", base, base + g3_len - 1)
        n_g3 += 1
        base += g3_len
        if base >= end_limit:
            break

        # 第一段 g1
        g1_len = predict_length("g1", base, stat, total_length,
                                length_models.get("g1"), fb_g1, min_val=1)
        if g1_len > max_g1:
            g1_len = max_g1
            cut_g1 += 1
        if base + g1_len - 1 > end_limit:
            break
        write_seg("g1", base, base + g1_len - 1)
        n_g1 += 1
        base += g1_len
        gene_blocks += 1  # 进入一个基因块

        # 该块内 (g2,g1) 对数
        pairs = predict_count(base, stat, total_length, count_model, fallback_count)
        if pairs > max_pairs:
            pairs = max_pairs
            cut_pairs += 1

        for _ in range(pairs):
            if base >= end_limit:
                break
            # g2
            g2_len = predict_length("g2", base, stat, total_length,
                                    length_models.get("g2"), fb_g2, min_val=10)
            if g2_len > max_g2:
                g2_len = max_g2
                cut_g2 += 1
            if base + g2_len - 1 > end_limit:
                break
            write_seg("g2", base, base + g2_len - 1)
            n_g2 += 1
            base += g2_len
            if base >= end_limit:
                break
            # 后续 g1
            g1_len = predict_length("g1", base, stat, total_length,
                                    length_models.get("g1"), fb_g1, min_val=1)
            if g1_len > max_g1:
                g1_len = max_g1
                cut_g1 += 1
            if base + g1_len - 1 > end_limit:
                break
            write_seg("g1", base, base + g1_len - 1)
            n_g1 += 1
            base += g1_len

    # 尾部若有余量再补一个 g3
    if base < end_limit:
        g3_len = predict_length("g3", base, stat, total_length,
                                length_models.get("g3"), fb_g3, min_val=g3min)
        if base + g3_len - 1 > end_limit:
            g3_len = end_limit - base + 1
        write_seg("g3", base, base + g3_len - 1)
        n_g3 += 1

    return {
        "n_g1": n_g1, "n_g2": n_g2, "n_g3": n_g3,
        "cut_g1": cut_g1, "cut_g2": cut_g2, "cut_pairs": cut_pairs,
        "gene_blocks": gene_blocks
    }

# -----------------------------
# 主程序
# -----------------------------
def main():
    ap = argparse.ArgumentParser(
        description="学习已有 SLiM 片段（-learn），并流式生成新的 SLiM 片段。"
    )
    ap.add_argument("-learn", action="store_true", help="启用学习模式")
    ap.add_argument("-i", "--input", required=True, help="输入 SLiM 片段文件")
    ap.add_argument("-l", "--length", type=int, required=True, help="染色体总长度（终止坐标）")
    ap.add_argument("-stat", type=int, default=1, help="起始坐标（默认 1）")
    ap.add_argument("-o", "--output", required=True, help="输出 SLiM 片段文件")
    ap.add_argument("-inv", required=True, help="学习区间，格式：start-end（如 100000000-200000000）")
    ap.add_argument("-ep", type=int, default=20, help="训练迭代次数（默认 20）")
    ap.add_argument("-g3min", type=int, default=100, help="g3 回退最小长度（默认 100）")
    ap.add_argument("-g3max", type=int, default=5000, help="g3 回退最大长度（默认 5000）")
    ap.add_argument("--g1-mean", type=float, help="g1 回退：原始尺度均值 E[X]")
    ap.add_argument("--g1-sigma", type=float, help="g1 回退：log 空间 σ")
    ap.add_argument("--g2-mean", type=float, help="g2 回退：原始尺度均值 E[X]")
    ap.add_argument("--g2-sigma", type=float, help="g2 回退：log 空间 σ")
    ap.add_argument("--max-g1", type=int, default=5000, help="g1 长度上限（默认 5000 bp）")
    ap.add_argument("--max-g2", type=int, default=200000, help="g2 长度上限（默认 200000 bp）")
    ap.add_argument("--max-pairs", type=int, default=100, help="每块最多 (g2,g1) 对数（默认 100）")
    ap.add_argument("--seed", type=int, help="随机种子（可选）")
    args = ap.parse_args()

    # 随机种子
    if args.seed is not None:
        np.random.seed(args.seed)
        random.seed(args.seed)
        tf.random.set_seed(args.seed)

    # 学习区间
    try:
        inv_start, inv_end = [int(x) for x in args.inv.split("-")]
    except Exception:
        sys.exit("学习区间 -inv 格式错误，应为 start-end（如 100000000-200000000）")
    inv_interval = (min(inv_start, inv_end), max(inv_start, inv_end))

    print("训练模型中……")
    length_models = train_length_models(args.input, inv_interval, args.stat, args.length, epochs=args.ep)
    count_model  = train_count_model(args.input, inv_interval, args.stat, args.length, epochs=args.ep)

    print("开始生成新的 SLiM 输出……")
    with open(args.output, "w") as fh:
        summary = deep_simulate_regions_stream(
            stat=args.stat,
            total_length=args.length,
            length_models=length_models,
            count_model=count_model,
            g3min=args.g3min,
            g3max=args.g3max,
            g1_mean=args.g1_mean,
            g1_sigma=args.g1_sigma,
            g2_mean=args.g2_mean,
            g2_sigma=args.g2_sigma,
            fh=fh,
            max_g1=args.max_g1,
            max_g2=args.max_g2,
            max_pairs=args.max_pairs
        )

    print(f"生成完成，新 SLiM 片段文件：{args.output}")
    print("—— 生成统计 ——")
    print(f"g1 段数: {summary['n_g1']}, 被上限截断: {summary['cut_g1']}")
    print(f"g2 段数: {summary['n_g2']}, 被上限截断: {summary['cut_g2']}")
    print(f"g3 段数: {summary['n_g3']}")
    print(f"基因块数: {summary['gene_blocks']}, pair 截断: {summary['cut_pairs']}")

if __name__ == "__main__":
    main()
