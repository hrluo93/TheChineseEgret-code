#!/usr/bin/env python3
"""
compute_g1_g2_g3_params.py

Estimate length distributions for exon (g1), intron (g2), and (optionally) intergenic (g3).
Optionally stratify by a categorical column (e.g., --type type).

Outputs:
  - For g1/g2 (log-normal):  --g1-mean <raw_mean> --g1-sigma <log_sigma>
                             --g2-mean <raw_mean> --g2-sigma <log_sigma>
  - For g3 (if provided):
      * Uniform (robust bounds from 10%–90% quantiles): g3min, g3max
      * Log-normal (optional): --g3-mean <raw_mean> --g3-sigma <log_sigma>

Usage:
  python compute_g1_g2_g3_params.py --exon EXON.csv --intron INTRON.csv [--inter INTER.csv]
         [--length-col LENGTH] [--trim 0.0] [--mle] [--type TYPECOL]

Notes
- If no `length` column exists, the script assumes cols 2/3 are start/end (0-based indexing for pandas),
  and computes length = end - start.
- By default uses sample std (ddof=1). Pass --mle to use MLE std (ddof=0).
- `--trim p` trims the top p fraction (e.g., 0.005) before computing raw mean for g1/g2.
- With --type TYPECOL, if TYPECOL exists in a table, stats are reported per category; otherwise falls back to overall.
"""

import argparse
import numpy as np
import pandas as pd

# ---------- I/O & preprocessing ----------

def read_table(path, length_col=None):
    """Read a table and return a DataFrame with a computed '_length' column."""
    sep = None
    if str(path).endswith((".tsv", ".bed")):
        sep = "\t"
    df = pd.read_csv(path, sep=sep, engine="python")
    cols = list(df.columns)

    if length_col and length_col in df.columns:
        lengths = df[length_col].astype(float)
    else:
        cand = [c for c in cols if c.lower() in ("length", "len", "size", "bp")]
        if cand:
            lengths = df[cand[0]].astype(float)
        else:
            if len(cols) < 3:
                raise ValueError(f"{path}: cannot find length column and fewer than 3 columns to derive from.")
            lengths = (df.iloc[:, 2] - df.iloc[:, 1]).astype(float)

    lengths = lengths.replace([np.inf, -np.inf], np.nan).dropna()
    lengths = lengths[lengths > 0]
    df = df.loc[lengths.index].copy()
    df["_length"] = lengths.astype(float)
    return df

# ---------- statistics & reporting ----------

def summarize_lengths(lengths, trim=0.0, ddof=1):
    n_all = len(lengths)
    raw_mean = float(lengths.mean())
    raw_median = float(lengths.median())
    raw_sd = float(lengths.std(ddof=ddof))
    raw_min = float(lengths.min())
    raw_max = float(lengths.max())

    # optional trimming on raw mean (for g1/g2; g3 uses untrimmed by default)
    if 0.0 < trim < 0.5:
        q = lengths.quantile(1.0 - trim)
        trimmed = lengths[lengths <= q]
        raw_mean_trim = float(trimmed.mean())
    else:
        trimmed = lengths
        raw_mean_trim = raw_mean

    # log-space stats (untrimmed by default)
    loglen = np.log(lengths)
    mu_hat = float(loglen.mean())
    sigma_hat = float(loglen.std(ddof=ddof))
    mean_check = float(np.exp(mu_hat + sigma_hat**2 / 2.0))

    # useful quantiles
    q05 = float(lengths.quantile(0.05))
    q10 = float(lengths.quantile(0.10))
    q90 = float(lengths.quantile(0.90))
    q95 = float(lengths.quantile(0.95))
    q99 = float(lengths.quantile(0.99))

    return {
        "n": n_all,
        "raw_mean": raw_mean,
        "raw_mean_trim": raw_mean_trim,
        "raw_median": raw_median,
        "raw_sd": raw_sd,
        "raw_min": raw_min,
        "raw_max": raw_max,
        "mu_hat": mu_hat,
        "sigma_hat": sigma_hat,
        "mean_check_from_mu_sigma": mean_check,
        "q05": q05, "q10": q10, "q90": q90, "q95": q95, "q99": q99,
        "trimmed_series": trimmed,  # kept for optional trimmed g3 estimates
    }

def pretty(name, s, prefix=None):
    title = name if prefix is None else f"{name} [{prefix}]"
    print(f"\n=== {title} ===")
    for k, v in s.items():
        if k == "trimmed_series":
            continue
        if isinstance(v, (float, np.floating)):
            print(f"{k:>26s}: {v:.6f}")
        else:
            print(f"{k:>26s}: {v}")

def recommend_args(tag, stats, use_trim=False, prefix=None):
    mean_used = stats["raw_mean_trim"] if use_trim else stats["raw_mean"]
    sigma = stats["sigma_hat"]
    mu_from_mean = np.log(mean_used) - 0.5 * sigma**2
    label = tag if prefix is None else f"{tag}[{prefix}]"
    print(f"\n-> Recommended arguments for {label} (log-normal):")
    print(f"   --{tag}-mean {mean_used:.6f}   --{tag}-sigma {sigma:.6f}")
    print(f"   (If sampling directly with μ/σ: μ = ln(mean) - 0.5·σ^2 = {mu_from_mean:.6f})")

def recommend_uniform_for_g3(stats, prefix=None):
    a = int(round(stats["q10"]))
    b = int(round(stats["q90"]))
    label = "g3" if prefix is None else f"g3[{prefix}]"
    print(f"\n-> Recommended bounds for {label} (Uniform):")
    print(f"   g3min = {a}   g3max = {b}   # based on 10%–90% quantiles")
    print("   (Also shown for reference: min/max and 5%/95% quantiles)")
    print(f"   raw min/max = {int(stats['raw_min'])} / {int(stats['raw_max'])}")
    print(f"   5% / 95%    = {int(round(stats['q05']))} / {int(round(stats['q95']))}")

def recommend_lognormal_for_g3(stats, prefix=None):
    mean_used = stats["raw_mean"]
    sigma = stats["sigma_hat"]
    mu_from_mean = np.log(mean_used) - 0.5 * sigma**2
    label = "g3" if prefix is None else f"g3[{prefix}]"
    print(f"\n-> Optional arguments for {label} (log-normal):")
    print(f"   --g3-mean {mean_used:.6f}   --g3-sigma {sigma:.6f}")
    print(f"   (μ from mean = {mu_from_mean:.6f}; caution: intergenic lengths are often heavy-tailed)")

def g3_trimmed_suggestions(lengths, ddof=1, prefix=None):
    for frac in (0.005, 0.01):  # 0.5% and 1% trims
        q = lengths.quantile(1 - frac)
        trimmed = lengths[lengths <= q]
        loglen_t = np.log(trimmed)
        mu_t = float(loglen_t.mean())
        sigma_t = float(loglen_t.std(ddof=ddof))
        mean_t = float(trimmed.mean())
        mu_from_mean_t = np.log(mean_t) - 0.5 * sigma_t**2
        label = "" if prefix is None else f" [{prefix}]"
        print(f"\n-> Optional (trim {frac*100:.1f}%) for g3{label} (log-normal, more robust):")
        print(f"   --g3-mean {mean_t:.6f}   --g3-sigma {sigma_t:.6f}")
        print(f"   (μ from mean = {mu_from_mean_t:.6f})")

# ---------- main ----------

def main():
    ap = argparse.ArgumentParser(
        description="Estimate log-normal parameters for g1/g2 and Uniform/Log-normal options for g3; optionally stratify by a categorical column."
    )
    ap.add_argument("--exon",   required=True, help="Exon (g1) file: CSV/TSV/BED-like.")
    ap.add_argument("--intron", required=True, help="Intron (g2) file: CSV/TSV/BED-like.")
    ap.add_argument("--inter",  help="Intergenic (g3) file: CSV/TSV/BED-like. Optional.")
    ap.add_argument("--length-col", help="Name of length column (if present).")
    ap.add_argument("--trim", type=float, default=0.0,
                    help="Trim top fraction before raw mean for g1/g2 (e.g., 0.005).")
    ap.add_argument("--mle", action="store_true",
                    help="Use MLE std (ddof=0) instead of sample std (ddof=1).")
    ap.add_argument("--type", dest="type_col",
                    help="Column name to stratify by (e.g., 'type'). If provided and present, compute per category.")
    args = ap.parse_args()

    ddof = 0 if args.mle else 1

    # ----- g1 (exon) -----
    exon_df = read_table(args.exon, args.length_col)
    if args.type_col and args.type_col in exon_df.columns:
        for cat, sub in exon_df.groupby(args.type_col, dropna=True):
            stats = summarize_lengths(sub["_length"], trim=args.trim, ddof=ddof)
            pretty("g1 (exon)", stats, prefix=str(cat))
            recommend_args("g1", stats, use_trim=(args.trim > 0), prefix=str(cat))
    else:
        stats = summarize_lengths(exon_df["_length"], trim=args.trim, ddof=ddof)
        pretty("g1 (exon)", stats)
        recommend_args("g1", stats, use_trim=(args.trim > 0))

    # ----- g2 (intron) -----
    intron_df = read_table(args.intron, args.length_col)
    if args.type_col and args.type_col in intron_df.columns:
        for cat, sub in intron_df.groupby(args.type_col, dropna=True):
            stats = summarize_lengths(sub["_length"], trim=args.trim, ddof=ddof)
            pretty("g2 (intron)", stats, prefix=str(cat))
            recommend_args("g2", stats, use_trim=(args.trim > 0), prefix=str(cat))
    else:
        stats = summarize_lengths(intron_df["_length"], trim=args.trim, ddof=ddof)
        pretty("g2 (intron)", stats)
        recommend_args("g2", stats, use_trim=(args.trim > 0))

    # ----- g3 (intergenic, optional) -----
    if args.inter:
        inter_df = read_table(args.inter, args.length_col)
        if args.type_col and args.type_col in inter_df.columns:
            for cat, sub in inter_df.groupby(args.type_col, dropna=True):
                stats = summarize_lengths(sub["_length"], trim=0.0, ddof=ddof)  # g3 默认不截尾
                pretty("g3 (intergenic)", stats, prefix=str(cat))
                recommend_uniform_for_g3(stats, prefix=str(cat))
                recommend_lognormal_for_g3(stats, prefix=str(cat))
                g3_trimmed_suggestions(sub["_length"], ddof=ddof, prefix=str(cat))
        else:
            stats = summarize_lengths(inter_df["_length"], trim=0.0, ddof=ddof)
            pretty("g3 (intergenic)", stats)
            recommend_uniform_for_g3(stats)
            recommend_lognormal_for_g3(stats)
            g3_trimmed_suggestions(inter_df["_length"], ddof=ddof)
    else:
        print("\n[Info] --inter not provided; skipping g3 (intergenic) estimation.")

if __name__ == "__main__":
    main()
