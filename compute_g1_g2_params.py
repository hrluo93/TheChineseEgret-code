#!/usr/bin/env python3
"""
compute_g1_g2_params.py

Estimate log-normal parameters for exon (g1) and intron (g2) lengths.
Outputs recommended arguments for your generator:
  --g1-mean <raw_mean> --g1-sigma <log_sigma>
  --g2-mean <raw_mean> --g2-sigma <log_sigma>

Usage:
  python compute_g1_g2_params.py --exon EXON.csv [--intron INTRON.csv]
         [--length-col LENGTH] [--trim 0.0] [--mle]

Notes
- If no `length` column exists, the script assumes cols 2/3 are start/end (0-based indexing for pandas),
  and computes length = end - start.
- By default uses sample std (ddof=1). Pass --mle to use MLE std (ddof=0).
- `--trim p` trims the top p fraction (e.g., 0.005) before computing raw mean to reduce long-tail influence.
"""

import argparse
import numpy as np
import pandas as pd
from pathlib import Path

def read_lengths(path, length_col=None):
    # Auto-detect separator
    sep = None
    if str(path).endswith((".tsv", ".bed")):
        sep = "\t"
    df = pd.read_csv(path, sep=sep, engine="python")
    cols = list(df.columns)
    if length_col and length_col in df.columns:
        lengths = df[length_col].astype(float)
    else:
        # Prefer explicit 'length' if present
        cand = [c for c in cols if c.lower() in ("length","len","size","bp")]
        if cand:
            lengths = df[cand[0]].astype(float)
        else:
            # Fall back to end - start using col2, col3
            if len(cols) < 3:
                raise ValueError(f"{path}: cannot find length column and fewer than 3 columns to derive from.")
            lengths = (df.iloc[:,2] - df.iloc[:,1]).astype(float)
    lengths = lengths.replace([np.inf, -np.inf], np.nan).dropna()
    # Keep positive lengths
    lengths = lengths[lengths > 0]
    return lengths

def summarize_lengths(lengths, trim=0.0, ddof=1):
    n_all = len(lengths)
    # raw stats
    raw_mean = float(lengths.mean())
    raw_median = float(lengths.median())
    raw_sd = float(lengths.std(ddof=ddof))

    # optional trimming on raw mean
    if 0.0 < trim < 0.5:
        q = lengths.quantile(1.0 - trim)
        trimmed = lengths[lengths <= q]
        raw_mean_trim = float(trimmed.mean())
    else:
        raw_mean_trim = raw_mean

    # log-space stats
    loglen = np.log(lengths)
    mu_hat = float(loglen.mean())
    sigma_hat = float(loglen.std(ddof=ddof))
    mean_check = float(np.exp(mu_hat + sigma_hat**2 / 2.0))

    return {
        "n": n_all,
        "raw_mean": raw_mean,
        "raw_mean_trim": raw_mean_trim,
        "raw_median": raw_median,
        "raw_sd": raw_sd,
        "mu_hat": mu_hat,
        "sigma_hat": sigma_hat,
        "mean_check_from_mu_sigma": mean_check
    }

def pretty(name, s):
    print(f"\n=== {name} ===")
    for k,v in s.items():
        print(f"{k:>26s}: {v:.6f}" if isinstance(v,(float,np.floating)) else f"{k:>26s}: {v}")

def recommend_args(tag, stats, use_trim=False):
    mean_used = stats["raw_mean_trim"] if use_trim else stats["raw_mean"]
    sigma = stats["sigma_hat"]
    mu_from_mean = np.log(mean_used) - 0.5 * sigma**2
    print(f"\n-> Recommended arguments for {tag}:")
    print(f"   --{tag}-mean {mean_used:.6f}   --{tag}-sigma {sigma:.6f}")
    print(f"   (If sampling directly with μ/σ: μ = ln(mean) - 0.5·σ^2 = {mu_from_mean:.6f})")

def main():
    ap = argparse.ArgumentParser(description="Estimate log-normal parameters for g1/g2.")
    ap.add_argument("--exon", required=True, help="Exon (g1) file: CSV/TSV/BED-like.")
    ap.add_argument("--intron", help="Intron (g2) file: CSV/TSV/BED-like.")
    ap.add_argument("--length-col", help="Name of length column (if present).")
    ap.add_argument("--trim", type=float, default=0.0, help="Trim top fraction before raw mean (e.g., 0.005).")
    ap.add_argument("--mle", action="store_true", help="Use MLE std (ddof=0) instead of sample std (ddof=1).")
    args = ap.parse_args()

    ddof = 0 if args.mle else 1

    # g1
    g1_len = read_lengths(args.exon, args.length_col)
    g1_stats = summarize_lengths(g1_len, trim=args.trim, ddof=ddof)
    pretty("g1 (exon)", g1_stats)
    recommend_args("g1", g1_stats, use_trim=(args.trim>0))

    # g2
    if args.intron:
        g2_len = read_lengths(args.intron, args.length_col)
        g2_stats = summarize_lengths(g2_len, trim=args.trim, ddof=ddof)
        pretty("g2 (intron)", g2_stats)
        recommend_args("g2", g2_stats, use_trim=(args.trim>0))

if __name__ == "__main__":
    main()
  
