#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compute region-length statistics and candidate distribution fits
for exon (g1), intron (g2), and intergenic (g3) files.

Key points
- Treat inputs as 1-based CLOSED intervals: length = end - start + 1
  (we ignore any existing "length" column if present).
- Optional category column (e.g., chromosome type) to report per-category fits.
- Fit candidate distributions (lognormal, gamma, exponential, normal, uniform)
  by MLE and compare with AIC.
- If best is lognormal, also print SLiM-friendly parameters:
    raw_mean (E[X]), sigma_hat (log-space), and μ = ln(mean) - 0.5·σ²

Usage
  python compute_region_params_1based.py \
      --exon   exon.bed \
      --intron intron.bed \
      --inter  inter.bed \
      --cat-col type \
      [--mle] [--trim 0.0]

Notes
- Separator auto-detected: if filename ends with .tsv/.bed => tab; else pandas auto.
- By default, log-space std uses sample std (ddof=1). Pass --mle to use ddof=0.
- --trim trims the top fraction before reporting *raw mean* (does not affect MLE fits).
"""

import argparse
import numpy as np
import pandas as pd
from scipy import stats

# ----------------------
# helpers
# ----------------------
def autodetect_sep(path_str: str):
    lower = path_str.lower()
    if lower.endswith(".tsv") or lower.endswith(".bed"):
        return "\t"
    return None  # let pandas infer

def read_table(path: str, cat_col: str | None):
    sep = autodetect_sep(path)
    df = pd.read_csv(path, sep=sep, engine="python")
    # normalize columns; we will compute length from 1-based CLOSED intervals
    cols = list(df.columns)
    if len(cols) < 3:
        raise ValueError(f"{path}: need at least 3 columns (chrom, start, end).")

    # start/end: use 2nd and 3rd columns regardless of their names
    start = pd.to_numeric(df.iloc[:, 1], errors="coerce")
    end   = pd.to_numeric(df.iloc[:, 2], errors="coerce")
    length = (end - start + 1).astype(float)

    out = pd.DataFrame({
        "length": length
    })

    if cat_col and (cat_col in df.columns):
        out["category"] = df[cat_col].astype(str)
    else:
        out["category"] = "ALL"

    # keep positive finite
    out = out.replace([np.inf, -np.inf], np.nan).dropna(subset=["length"])
    out = out[out["length"] > 0]
    return out

def summarize_lengths(arr: np.ndarray, trim: float = 0.0, ddof: int = 1):
    arr = np.asarray(arr, dtype=float)
    n = arr.size
    q1 = float(np.quantile(arr, 0.25))
    q3 = float(np.quantile(arr, 0.75))
    mean_raw = float(np.mean(arr))
    sd_raw = float(np.std(arr, ddof=ddof))
    median = float(np.median(arr))
    vmin = float(np.min(arr))
    vmax = float(np.max(arr))

    if 0.0 < trim < 0.5:
        thr = np.quantile(arr, 1.0 - trim)
        mean_trim = float(np.mean(arr[arr <= thr]))
    else:
        mean_trim = mean_raw

    # log-space estimators (sample by default)
    loglen = np.log(arr)
    mu_hat = float(np.mean(loglen))
    sigma_hat = float(np.std(loglen, ddof=ddof))
    mean_check = float(np.exp(mu_hat + 0.5 * sigma_hat**2))

    return {
        "n": n,
        "mean": mean_raw,
        "mean_trim": mean_trim,
        "median": median,
        "sd": sd_raw,
        "q1": q1,
        "q3": q3,
        "min": vmin,
        "max": vmax,
        "mu_hat": mu_hat,          # log-space mean (empirical)
        "sigma_hat": sigma_hat,    # log-space std (empirical)
        "mean_check": mean_check   # exp(mu + sigma^2 / 2)
    }

def aic_from_loglik(loglik: float, k_params: int) -> float:
    return 2 * k_params - 2 * loglik

def fit_and_aic(arr: np.ndarray):
    """
    Fit candidate distributions by MLE and compute AIC.
    Returns dict keyed by dist name with {params, aic, loglik, k}.
    - lognorm: (s, loc=0, scale)  => k = 2 (s, scale)  [loc fixed]
    - gamma  : (a, loc=0, scale)  => k = 2 (a, scale)  [loc fixed]
    - expon  : (loc=0, scale)     => k = 1 (scale)     [loc fixed]
    - norm   : (loc, scale)       => k = 2
    - uniform: (loc, scale)       => k = 2
    """
    out = {}
    x = np.asarray(arr, dtype=float)

    # Lognormal
    try:
        s, loc, scale = stats.lognorm.fit(x, floc=0)
        dist = stats.lognorm(s, loc=loc, scale=scale)
        ll = float(np.sum(dist.logpdf(x)))
        out["lognorm"] = {"params": (s, loc, scale), "loglik": ll, "k": 2,
                          "aic": aic_from_loglik(ll, 2)}
    except Exception:
        pass

    # Gamma
    try:
        a, loc, scale = stats.gamma.fit(x, floc=0)
        dist = stats.gamma(a, loc=loc, scale=scale)
        ll = float(np.sum(dist.logpdf(x)))
        out["gamma"] = {"params": (a, loc, scale), "loglik": ll, "k": 2,
                        "aic": aic_from_loglik(ll, 2)}
    except Exception:
        pass

    # Exponential
    try:
        loc, scale = stats.expon.fit(x, floc=0)
        dist = stats.expon(loc=loc, scale=scale)
        ll = float(np.sum(dist.logpdf(x)))
        out["expon"] = {"params": (loc, scale), "loglik": ll, "k": 1,
                        "aic": aic_from_loglik(ll, 1)}
    except Exception:
        pass

    # Normal (not ideal for strictly-positive heavy-tailed data, but for completeness)
    try:
        loc, scale = stats.norm.fit(x)
        dist = stats.norm(loc=loc, scale=scale)
        ll = float(np.sum(dist.logpdf(x)))
        out["norm"] = {"params": (loc, scale), "loglik": ll, "k": 2,
                       "aic": aic_from_loglik(ll, 2)}
    except Exception:
        pass

    # Uniform
    try:
        loc, scale = stats.uniform.fit(x)
        dist = stats.uniform(loc=loc, scale=scale)
        ll = float(np.sum(dist.logpdf(x)))
        out["uniform"] = {"params": (loc, scale), "loglik": ll, "k": 2,
                          "aic": aic_from_loglik(ll, 2)}
    except Exception:
        pass

    return out

def print_block(title: str):
    print("\n" + title)
    print("-" * len(title))

def report_one(name: str, df: pd.DataFrame, trim: float, ddof: int):
    if df.empty:
        print(f"\n===== {name} =====\n(No data)")
        return

    print(f"\n===== {name} =====")
    for cat, sub in df.groupby("category", sort=False):
        arr = sub["length"].to_numpy()
        stats_basic = summarize_lengths(arr, trim=trim, ddof=ddof)

        print(f"\n--- Category: {cat} ---")
        print(f"n={stats_basic['n']:,}, mean={stats_basic['mean']:.6f}, "
              f"median={stats_basic['median']:.6f}, sd={stats_basic['sd']:.6f}, "
              f"Q1={stats_basic['q1']:.6f}, Q3={stats_basic['q3']:.6f}, "
              f"min={stats_basic['min']:.6f}, max={stats_basic['max']:.6f}")

        # Fit distributions
        fits = fit_and_aic(arr)
        if not fits:
            print("No fit produced (check data).")
            continue

        # choose best by AIC
        best_name = min(fits.keys(), key=lambda k: fits[k]["aic"])
        best = fits[best_name]
        print(f"Best by AIC: {best_name}  (AIC={best['aic']:.2f})")

        # pretty print parameters depending on dist
        if best_name == "lognorm":
            s, loc, scale = best["params"]
            # For lognorm in SciPy: shape=s(=σ), loc, scale=exp(μ)
            sigma_fit = s
            mu_fit = np.log(scale)
            ex_mean_fit = float(np.exp(mu_fit + 0.5 * sigma_fit**2))
            print(f"  lognorm params (shape=σ, loc={loc:.0f}, scale=e^μ): σ={sigma_fit:.6f}, μ={mu_fit:.6f}, E[X]={ex_mean_fit:.6f}")

            # Also give log-space empirical sigma for SLiM fallback
            sigma_emp = stats_basic["sigma_hat"]
            mean_raw = stats_basic["mean_trim"] if trim > 0 else stats_basic["mean"]
            mu_from_mean = np.log(mean_raw) - 0.5 * sigma_emp**2
            check_back = float(np.exp(mu_from_mean + 0.5 * sigma_emp**2))
            print("\n# Lognormal-friendly recommendation (for SLiM fallback):")
            print(f"   raw_mean (E[X])        : {mean_raw:.6f}")
            print(f"   sigma_hat (log-space σ): {sigma_emp:.6f}")
            print(f"   μ = ln(mean) - 0.5·σ²  : {mu_from_mean:.6f}")
            print(f"   (check E[X] from μ/σ   : {check_back:.6f})")

        elif best_name == "gamma":
            a, loc, scale = best["params"]
            print(f"  gamma params (a, loc={loc:.0f}, scale): a={a:.6f}, scale={scale:.6f}")
        elif best_name == "expon":
            loc, scale = best["params"]
            print(f"  expon params (loc={loc:.0f}, scale): scale={scale:.6f}")
        elif best_name == "norm":
            loc, scale = best["params"]
            print(f"  norm params (mean, sd): mean={loc:.6f}, sd={scale:.6f}")
        elif best_name == "uniform":
            loc, scale = best["params"]
            print(f"  uniform params (loc, scale): loc={loc:.6f}, scale={scale:.6f} "
                  f"[support=({loc:.2f}, {loc+scale:.2f})]")

# ----------------------
# main
# ----------------------
def main():
    ap = argparse.ArgumentParser(description="Estimate region-length distributions (1-based CLOSED intervals).")
    ap.add_argument("--exon",   help="Exon file (CSV/TSV/BED-like).")
    ap.add_argument("--intron", help="Intron file (CSV/TSV/BED-like).")
    ap.add_argument("--inter",  help="Intergenic file (CSV/TSV/BED-like).")
    ap.add_argument("--cat-col", help="Optional category column (e.g., chromosome type).")
    ap.add_argument("--trim", type=float, default=0.0, help="Trim top fraction before raw mean (e.g., 0.005).")
    ap.add_argument("--mle", action="store_true", help="Use MLE std (ddof=0) instead of sample std (ddof=1).")
    args = ap.parse_args()

    ddof = 0 if args.mle else 1
    any_file = False

    if args.exon:
        any_file = True
        df = read_table(args.exon, args.cat_col)
        report_one("g1 (exon)", df, trim=args.trim, ddof=ddof)

    if args.intron:
        any_file = True
        df = read_table(args.intron, args.cat_col)
        report_one("g2 (intron)", df, trim=args.trim, ddof=ddof)

    if args.inter:
        any_file = True
        df = read_table(args.inter, args.cat_col)
        report_one("g3 (intergenic)", df, trim=args.trim, ddof=ddof)

    if not any_file:
        raise SystemExit("Please provide at least one of --exon / --intron / --inter.")

if __name__ == "__main__":
    main()