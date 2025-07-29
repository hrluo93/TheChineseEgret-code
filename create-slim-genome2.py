#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
学习模式脚本（-learn）：
从已有 SLiM initializeGenomicElement(...) 文件中，在指定学习区间内拟合
g1/g2/g3 的长度与 (g2,g1) 对的数量（基因块大小）与位置的关系（简单 MLP 回归），
再按 g3→g1→(g2,g1)*→g3 生成新的片段序列。

新增：
- --g1-mean/--g1-sigma/--g2-mean/--g2-sigma：当模型不可用时，作为 g1/g2 的
  对数正态 fallback 参数（原始尺度均值 + log 空间 σ）。未提供时退回“手册示例”：
  g1 ~ LogNormal(ln(50), ln(2)), g2 ~ LogNormal(ln(100), ln(1.5)).
- --seed：设定随机种子，保证复现。
"""
import argparse
import random
import sys
import re
import numpy as np
import tensorflow as tf
from tensorflow import keras
from keras import layers, models

##############################
# 训练数据提取
##############################
def parse_slim_line(line):
    """
    从格式类似：
    initializeGenomicElement(g3, 188618852, 188693529);
    中提取 (feature, start, end)
    """
    pattern = r'initializeGenomicElement\(\s*([^,]+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\);'
    m = re.match(pattern, line)
    if m:
        feature = m.group(1).strip()
        start = int(m.group(2))
        end = int(m.group(3))
        return (feature, start, end)
    return None

def read_length_data(input_file, inv_interval, stat, total_length):
    """
    从 slim 文件中提取长度数据：对于每条记录（g1, g2, g3），
    如果该记录的中点在 inv_interval 内，则记录 (x, length)
    其中 x = (mid - stat)/total_length 为归一化位置，length = end - start + 1
    返回字典 {feature: [(x, length), ...]}
    """
    data = {"g1": [], "g2": [], "g3": []}
    inv_start, inv_end = inv_interval
    with open(input_file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            seg = parse_slim_line(line)
            if seg:
                feature, s, e = seg
                mid = (s + e) // 2
                if inv_start <= mid <= inv_end:
                    x = (mid - stat) / total_length
                    length = e - s + 1
                    if feature in data:
                        data[feature].append((x, length))
    return data

def read_count_data(input_file, inv_interval, stat, total_length):
    """
    将 slim 文件按基因块分组。假设格式为：
       g3 (noncoding boundary)
       g1 (first exon)
       then 反复： g2, g1 的对
       最后一个 g3 作为边界。
    对于每个基因块（即从一个 g3 到下一个 g3），取块内（除边界）的记录，
    统计 (g2, g1) 对的数量。并以该块代表位置（第一条 g1 的中点）归一化作为输入 x。
    返回列表 [(x, count), ...]
    """
    blocks = []
    with open(input_file) as f:
        lines = [l.strip() for l in f if l.strip()]
    # 解析所有行
    segs = [parse_slim_line(line) for line in lines if parse_slim_line(line) is not None]
    if not segs:
        return blocks
    # 寻找 g3 的索引
    indices = [i for i, seg in enumerate(segs) if seg[0] == "g3"]
    # 每个块：从 indices[i] 到 indices[i+1]（不包含后者）
    for i in range(len(indices) - 1):
        block = segs[indices[i]:indices[i + 1]]
        # 如果块内至少有3条记录（g3, g1, g3），则有基因信息
        if len(block) < 3:
            continue
        # 假设第一条为 g3，最后一条为 g3，块内剩余为基因信息
        gene_block = block[1:-1]
        # 验证 gene_block：第一条必须为 g1，其余应按 (g2, g1) 成对出现
        if gene_block[0][0] != "g1":
            continue
        pair_count = 0
        # 从第二个记录开始，每两个记录应为 (g2, g1)
        j = 1
        while j < len(gene_block) - 1:
            if gene_block[j][0] == "g2" and gene_block[j + 1][0] == "g1":
                pair_count += 1
                j += 2
            else:
                j += 1
        # 定义块的代表位置：使用第一条 g1 中点
        g1_seg = gene_block[0]
        g1_mid = (g1_seg[1] + g1_seg[2]) // 2
        if inv_interval[0] <= g1_mid <= inv_interval[1]:
            x = (g1_mid - stat) / total_length
            blocks.append((x, pair_count))
    return blocks

##############################
# 模型训练
##############################
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
    models_dict = {}
    for feature in ["g1", "g2", "g3"]:
        feat_data = data.get(feature, [])
        if len(feat_data) == 0:
            print(f"Warning: No training data for {feature} in learning interval!")
            models_dict[feature] = None
        else:
            arr = np.array(feat_data)  # shape (n,2)
            X = arr[:, 0].reshape(-1, 1).astype(np.float32)
            y = arr[:, 1].reshape(-1, 1).astype(np.float32)
            model = build_regression_model()
            model.fit(X, y, epochs=epochs, verbose=0)
            preds = model.predict(X, verbose=0)
            rmse = np.sqrt(np.mean((preds.flatten() - y.flatten()) ** 2))
            print(f"Feature {feature}: trained on {len(X)} samples, RMSE = {rmse:.2f}")
            models_dict[feature] = (model, rmse)
    return models_dict

def train_count_model(input_file, inv_interval, stat, total_length, epochs=20):
    blocks = read_count_data(input_file, inv_interval, stat, total_length)
    if len(blocks) == 0:
        print("Warning: No count training data in learning interval!")
        return None
    arr = np.array(blocks)  # shape (n,2)
    X = arr[:, 0].reshape(-1, 1).astype(np.float32)
    y = arr[:, 1].reshape(-1, 1).astype(np.float32)  # count
    model = build_regression_model()
    model.fit(X, y, epochs=epochs, verbose=0)
    preds = model.predict(X, verbose=0)
    rmse = np.sqrt(np.mean((preds.flatten() - y.flatten()) ** 2))
    print(f"Count model: trained on {len(X)} samples, RMSE = {rmse:.2f}")
    return (model, rmse)

def train_models(input_file, inv_interval, stat, total_length, epochs=20):
    length_models = train_length_models(input_file, inv_interval, stat, total_length, epochs)
    count_model = train_count_model(input_file, inv_interval, stat, total_length, epochs)
    return length_models, count_model

##############################
# 预测函数：长度与数量
##############################
def predict_length(feature, base, stat, total_length, model_info, fallback_value, min_val):
    """
    对于 feature，根据当前 base 计算归一化输入，
    若模型可用，则预测长度并加上噪声（标准差为模型 RMSE），否则使用 fallback_value()。
    返回至少 min_val。
    """
    if model_info is not None:
        model, rmse = model_info
        x = np.array([[(base - stat) / total_length]], dtype=np.float32)
        pred = float(model.predict(x, verbose=0)[0, 0])
        noise = np.random.normal(0, rmse)
        length = int(round(pred + noise))
    else:
        length = int(fallback_value())
    return max(length, min_val)

def predict_count(base, stat, total_length, model_info, fallback_fn):
    """
    预测 count（内含子–外显子对数量）的函数：
    若模型可用，则用回归模型预测，并取 round；否则使用 fallback_fn 计算数量。
    """
    if model_info is not None:
        model, rmse = model_info
        x = np.array([[(base - stat) / total_length]], dtype=np.float32)
        pred = float(model.predict(x, verbose=0)[0, 0])
        noise = np.random.normal(0, rmse)
        cnt = int(round(pred + noise))
    else:
        cnt = int(fallback_fn())
    return max(cnt, 0)

def fallback_count():
    """
    fallback: 使用几何分布模拟 count，即不断累计直到失败 (概率 0.8 继续生成)
    """
    cnt = 0
    while random.random() < 0.8:
        cnt += 1
    return cnt

##############################
# 模拟生成：使用深度学习预测长度与数量
##############################
def deep_simulate_regions(stat, total_length, length_models, count_model, g3min, g3max,
                          g1_mean=None, g1_sigma=None, g2_mean=None, g2_sigma=None):
    """
    深度学习模式下模拟生成 slim 调用：
      - g3（非编码）：采用 predict_length，fallback 为均匀采样 [g3min, g3max]
      - g1（外显子）：采用 predict_length，fallback：
            若提供 g1_mean/g1_sigma 则 ~ LogNormal(mu_ex, g1_sigma)，
            其中 mu_ex = ln(g1_mean) - 0.5*g1_sigma^2；否则用“手册默认” ln(50), ln(2)。
      - g2（内含子）：同上；若提供 g2_mean/g2_sigma 则按参数采样，否则用“手册默认” ln(100), ln(1.5)。
      - gene 块中数量由 predict_count (fallback 使用 fallback_count) 预测
    模拟流程：
      while(base < total_length):
          生成一段 g3
          如果 base < total_length:
              生成第一段 g1
              预测 count
              for i in range(count):
                  生成 g2
                  生成后续 g1
      最后若剩余生成一段 g3
    返回模拟区段列表，每个为 (feature, start, end)
    """
    sim_segments = []
    base = stat

    # fallback for g3
    fb_g3 = lambda: random.randint(g3min, g3max)

    # fallback for g1
    if (g1_mean is not None) and (g1_sigma is not None):
        mu_ex = np.log(float(g1_mean)) - 0.5 * (float(g1_sigma) ** 2)
        fb_g1 = lambda: int(np.random.lognormal(mu_ex, float(g1_sigma))) + 1
    else:
        # 手册默认
        fb_g1 = lambda: int(np.random.lognormal(np.log(50.0), np.log(2.0))) + 1

    # fallback for g2
    if (g2_mean is not None) and (g2_sigma is not None):
        mu_in = np.log(float(g2_mean)) - 0.5 * (float(g2_sigma) ** 2)
        fb_g2 = lambda: int(np.random.lognormal(mu_in, float(g2_sigma))) + 10
    else:
        # 手册默认
        fb_g2 = lambda: int(np.random.lognormal(np.log(100.0), np.log(1.5))) + 10

    while base < total_length:
        # 非编码区 (g3)
        g3_len = predict_length("g3", base, stat, total_length, length_models.get("g3"), fb_g3, g3min)
        if base + g3_len - 1 > total_length:
            g3_len = total_length - base + 1
        sim_segments.append(("g3", base, base + g3_len - 1))
        base += g3_len
        if base >= total_length:
            break

        # 第一个外显子 (g1)
        g1_len = predict_length("g1", base, stat, total_length, length_models.get("g1"), fb_g1, 1)
        if base + g1_len - 1 > total_length:
            break
        sim_segments.append(("g1", base, base + g1_len - 1))
        base += g1_len

        # 预测该 gene 块中 (g2, g1) 对的数量
        count = predict_count(base, stat, total_length, count_model, fallback_count)
        for _ in range(count):
            if base >= total_length:
                break
            g2_len = predict_length("g2", base, stat, total_length, length_models.get("g2"), fb_g2, 10)
            if base + g2_len - 1 > total_length:
                break
            sim_segments.append(("g2", base, base + g2_len - 1))
            base += g2_len
            if base >= total_length:
                break
            g1_len = predict_length("g1", base, stat, total_length, length_models.get("g1"), fb_g1, 1)
            if base + g1_len - 1 > total_length:
                break
            sim_segments.append(("g1", base, base + g1_len - 1))
            base += g1_len

    # 尾部补一段 g3
    if base < total_length:
        g3_len = predict_length("g3", base, stat, total_length, length_models.get("g3"), fb_g3, g3min)
        if base + g3_len - 1 > total_length:
            g3_len = total_length - base + 1
        sim_segments.append(("g3", base, base + g3_len - 1))

    return sim_segments

##############################
# 输出结果
##############################
def output_simulation(sim_segments, out_file):
    lines = []
    for seg in sim_segments:
        feat, s, e = seg
        lines.append(f"initializeGenomicElement({feat}, {s}, {e});")
    with open(out_file, "w") as f:
        f.write("\n".join(lines) + "\n")

##############################
# 主程序
##############################
def main():
    parser = argparse.ArgumentParser(
        description="学习模式脚本（-learn）：通过输入已有 slim 文件学习各区段长度和 g1/g2 数量分布，"
                    "并利用深度学习模型生成新的 slim 输出。\n"
                    "参数 -l 指定染色体总长度，-stat 指定起始坐标，-o 输出文件，-inv 指定学习区间（格式 start-end）。"
    )
    parser.add_argument("-learn", action="store_true", help="启用学习模式")
    parser.add_argument("-i", "--input", required=True, help="输入已有 slim 文件")
    parser.add_argument("-l", "--length", type=int, required=True, help="染色体模拟总长度")
    parser.add_argument("-stat", type=int, default=1, help="起始坐标，默认 1")
    parser.add_argument("-o", "--output", required=True, help="输出 slim 文件")
    parser.add_argument("-inv", required=True, help="学习区间，格式如：100000000-200000000")
    parser.add_argument("-ep", type=int, default=20, help="训练迭代次数，默认 20")
    parser.add_argument("-g3min", type=int, default=100, help="fallback 非编码区最小长度，默认 100")
    parser.add_argument("-g3max", type=int, default=5000, help="fallback 非编码区最大长度，默认 5000")
    # 新增：可选的 g1/g2 fallback 参数（原始均值 + log σ）
    parser.add_argument("--g1-mean", type=float, help="fallback g1 原始均值 E[X]（可选）")
    parser.add_argument("--g1-sigma", type=float, help="fallback g1 log 空间 σ（可选）")
    parser.add_argument("--g2-mean", type=float, help="fallback g2 原始均值 E[X]（可选）")
    parser.add_argument("--g2-sigma", type=float, help="fallback g2 log 空间 σ（可选）")
    # 新增：随机种子
    parser.add_argument("--seed", type=int, help="随机种子（可选）")

    args = parser.parse_args()

    # 设定随机种子，保证复现（可选）
    if args.seed is not None:
        np.random.seed(args.seed)
        random.seed(args.seed)
        tf.random.set_seed(args.seed)

    total_length = args.length
    stat = args.stat
    try:
        inv_start, inv_end = [int(x) for x in args.inv.split("-")]
    except Exception:
        sys.exit("学习区间 -inv 格式错误，应为: start-end, 例如 100000000-200000000")
    inv_interval = (inv_start, inv_end)

    print("训练模型中……")
    length_models, count_model = train_models(args.input, inv_interval, stat, total_length, epochs=args.ep)

    print("开始生成新的 slim 输出……")
    sim_segments = deep_simulate_regions(
        stat, total_length, length_models, count_model, args.g3min, args.g3max,
        g1_mean=args.g1_mean, g1_sigma=args.g1_sigma,
        g2_mean=args.g2_mean, g2_sigma=args.g2_sigma
    )
    output_simulation(sim_segments, args.output)
    print(f"生成完成，新 slim 输出文件为 {args.output}")

if __name__ == "__main__":
    main()
