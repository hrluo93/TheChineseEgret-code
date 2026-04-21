#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
easySFS 1D (.sfs / .obs) -> mushi folded-SFS to η(t)
输出：
  <out>.ksfs.tsv        # folded SFS (freq\tcount)
  <out>.eta.tsv         #  η(t)=2N(t)
  <out>.eta.svg/.pdf    # Output
  <out>.overlay.svg/pdf # （--grid）多组正则结果叠加图
"""

import argparse, os, re
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import mushi

# ---------- 解析 easySFS 文件 ----------
def parse_dadi_sfs(path):
    lines = []
    with open(path, 'r') as f:
        for ln in f:
            s = ln.strip()
            if s and not s.startswith('#'):
                lines.append(s)
    if len(lines) < 3:
        raise ValueError("Unexpected dadi .sfs format (need dims/data/mask).")
    header = lines[0].split()
    dims = [int(x) for x in header if x.isdigit()]
    D = dims[0]
    is_folded = any(tok.lower() == 'folded' for tok in header)
    data = np.array([float(x) for x in lines[1].split()], dtype=float)
    mask = np.array([int(x) for x in lines[2].split()], dtype=int)
    data = np.where(mask == 1, 0.0, data)
    if is_folded:
        # 从 singleton 开始（去掉 monomorphic 0）
        counts = data[1:].astype(float)
        freqs  = np.arange(1, len(counts) + 1)
    else:
        # 未折叠：折叠为 MAF（不含 0 和 n）
        n = D - 1
        kmax = n // 2
        folded = []
        for i in range(1, kmax + 1):
            folded.append(data[i] if i == n - i else data[i] + data[n - i])
        counts = np.array(folded, dtype=float)
        freqs  = np.arange(1, kmax + 1)
    return freqs, counts

def parse_fsc_obs(path):
    with open(path, 'r') as f:
        lines = [ln.strip() for ln in f if ln.strip()]
    nums = np.array([float(x) for x in re.findall(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', lines[1])], dtype=float)
    name = os.path.basename(path).lower()
    is_folded = ('mafpop' in name) or ('maf' in name)
    if is_folded:
        # 从 singleton 开始（去掉 monomorphic 0）
        counts = nums[1:]
        freqs  = np.arange(1, len(nums))
    else:
        n = len(nums) - 1
        kmax = n // 2
        folded = []
        for i in range(1, kmax + 1):
            folded.append(nums[i] if i == n - i else nums[i] + nums[n - i])
        counts = np.array(folded, dtype=float)
        freqs  = np.arange(1, kmax + 1)
    return freqs, counts

def write_ksfs_tsv(path, freqs, counts):
    with open(path, 'w') as w:
        w.write("freq\tcount\n")
        for f, c in zip(freqs, counts):
            w.write(f"{int(f)}\t{int(round(c))}\n")

# ---------- 解析 penalty ----------
def parse_one_penalty_set(s):
    """'0,300;1,12;2,1' -> [(0,300.0),(1,12.0),(2,1.0)]"""
    s = s.strip()
    if not s:
        return []
    out = []
    for seg in s.split(";"):
        seg = seg.strip()
        if not seg:
            continue
        a, b = seg.split(",")
        out.append((int(a), float(b)))
    return out

def penalty_label(pens):
    return "p_" + "_".join(f"{o}-{int(v) if v.is_integer() else v}" for o, v in pens)

# ---------- 主流程 ----------
def main():
    ap = argparse.ArgumentParser(description="easySFS -> mushi (folded) -> SVG/PDF（AI 可编辑）")
    ap.add_argument("-i", "--input", required=True, help=".sfs (dadi) 或 .obs (fastsimcoal2)")
    ap.add_argument("-o", "--outprefix", required=True)
    ap.add_argument("-n", "--bins", type=int, required=True,
                    help="手动指定折叠频谱保留的 bin 数（从 freq=1 起，保留前 n 行）")
    ap.add_argument("--mu", type=float, required=True, help="每碱基每代突变率 (e.g. 1.25e-8)")
    ap.add_argument("--L",  type=float, required=True, help="callable 基因组长度 (bp)")
    ap.add_argument("--penalty", default="0,300;1,12",
                    help='单组趋势罚项（如 "0,300;1,12"）')
    ap.add_argument("--grid", default=None,
                    help='多组趋势罚项，用 | 分隔：如 "0,300;1,12|0,200;1,25;2,1"')
    ap.add_argument("--ridge", type=float, default=1e-4, help="ridge_penalty（默认 1e-4）")
    ap.add_argument("--min-freq", type=int, default=1,
                    help="最小频率类（默认1；设2屏蔽singleton）")
    ap.add_argument("--max-iter", type=int, default=300)
    ap.add_argument("--format", choices=["svg","pdf"], default="svg", help="矢量输出格式（默认 svg）")
    ap.add_argument("--font", default=None, help="首选字体名（如 Noto Sans CJK SC；需系统已安装）")
    ap.add_argument("--transparent", action="store_true", help="透明背景导出")
    args = ap.parse_args()

    # 让文字在 SVG/PDF 中可编辑
    mpl.rcParams['text.usetex'] = False
    if args.format == "svg":
        mpl.rcParams['svg.fonttype'] = 'none'
    else:
        mpl.rcParams['pdf.fonttype'] = 42
        mpl.rcParams['ps.fonttype']  = 42
    if args.font:
        mpl.rcParams['font.family'] = 'sans-serif'
        mpl.rcParams['font.sans-serif'] = [args.font]

    # 读取 SFS
    ext = os.path.splitext(args.input)[1].lower()
    if ext == ".sfs":
        freqs, counts = parse_dadi_sfs(args.input)
    elif ext == ".obs":
        freqs, counts = parse_fsc_obs(args.input)
    else:
        try:
            freqs, counts = parse_dadi_sfs(args.input)
        except Exception:
            freqs, counts = parse_fsc_obs(args.input)
    if counts.sum() == 0:
        raise SystemExit("All counts are zero after parsing; check input file.")

    # 写 mushi 可读的两列表 TSV（保留表头）
    tsv = f"{args.outprefix}.ksfs.tsv"
    write_ksfs_tsv(tsv, freqs, counts)

    # 读回并裁剪（容错：有/无表头都能读）
    try:
        df = pd.read_csv(tsv, sep="\t", comment="#", engine="python")
        if not {"freq","count"}.issubset(df.columns):
            df = pd.read_csv(tsv, sep="\t", header=None,
                             names=["freq","count"], comment="#", engine="python")
    except Exception:
        df = pd.read_csv(tsv, sep="\t", header=None,
                         names=["freq","count"], comment="#", engine="python")

    # 频率筛选与截取
    df["freq"] = df["freq"].astype(int)
    df["count"] = df["count"].astype(float)
    df = df[df["freq"] >= args.min_freq].copy()
    if len(df) < args.bins:
        raise SystemExit(f"请求保留 n={args.bins} 个 bin，但可用 bin 仅 {len(df)}（受 --min-freq 影响）。")
    df = df.iloc[:args.bins].copy()

    X = df["count"].to_numpy()
    mu0 = args.mu * args.L

    # 解析 grid / 单组
    grid_sets = []
    if args.grid:
        for chunk in args.grid.split("|"):
            pens = parse_one_penalty_set(chunk)
            if pens:
                grid_sets.append(pens)
    else:
        grid_sets = [parse_one_penalty_set(args.penalty)]

    # 叠加图
    overlay = (len(grid_sets) > 1)
    if overlay:
        plt.figure()

    results = []
    for pens in grid_sets:
        label = penalty_label(pens)
        # 初始化新的 ksfs，避免状态干扰
        ks = mushi.kSFS(X=X)
        ks.infer_eta(mu0, *pens,
                     ridge_penalty=args.ridge,
                     max_iter=args.max_iter,
                     verbose=True, folded=True)

        # 导出每组的 eta.tsv 和 单独曲线
        tag = f"{args.outprefix}.{label}"
        with open(f"{tag}.eta.tsv", "w") as w:
            w.write("t_start\tt_end\teta_2N\n")
            for t0, t1, val in ks.eta.epochs():
                w.write(f"{t0}\t{t1}\t{val}\n")

        plt.figure()
        ks.eta.plot(label="mushi (folded)")
        plt.xscale("log"); plt.yscale("log")
        plt.xlabel("Generations ago"); plt.ylabel("η(t) = 2N(t)")
        plt.legend(); plt.tight_layout()
        plt.savefig(f"{tag}.eta.{args.format}", format=args.format,
                    bbox_inches="tight", transparent=args.transparent)

        if overlay:
            # 叠加到同一张图（当前激活的 overlay fig）
            ks.eta.plot(label=label)

        results.append(label)

    if overlay:
        plt.xscale("log"); plt.yscale("log")
        plt.xlabel("Generations ago"); plt.ylabel("η(t) = 2N(t)")
        plt.legend(); plt.tight_layout()
        plt.savefig(f"{args.outprefix}.overlay.{args.format}", format=args.format,
                    bbox_inches="tight", transparent=args.transparent)

    print("Ran penalty sets:", ", ".join(results))
    print("Wrote:", tsv)
    if overlay:
        print(f"Overlay figure: {args.outprefix}.overlay.{args.format}")

if __name__ == "__main__":
    main()
  
