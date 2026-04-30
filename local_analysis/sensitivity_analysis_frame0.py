#!/usr/bin/env python3
"""
sensitivity_analysis_frame0.py — LOCAL USE ONLY, do not push to GitHub.

Shows how ORF yield changes across a range of frame0_fraction thresholds
and min_psite_reads values. Purpose: justify the hard filter values used
in the paper (frame0 >= 0.55, min_psites >= 8).

Run this on your Stage 5 "with_periodicity" output (the full table before
the hard filter is applied) to see the sensitivity curve.

Usage:
    python sensitivity_analysis_frame0.py \
        --input   results/phase4/stage5_HeLa_M_rep1_with_periodicity.csv \
        --sample  HeLa_M_rep1 \
        --outdir  local_analysis/

Output:
    - sensitivity_HeLa_M_rep1.csv      full threshold sweep table
    - sensitivity_HeLa_M_rep1.pdf      4-panel figure
"""

import os
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--input",   required=True, help="Stage 5 with_periodicity CSV")
    p.add_argument("--sample",  required=True)
    p.add_argument("--outdir",  default="local_analysis")
    p.add_argument("--frame_range", nargs=2, type=float, default=[0.33, 0.85],
                   help="Min and max frame0_fraction to sweep")
    p.add_argument("--psite_range", nargs=2, type=int, default=[2, 20],
                   help="Min and max min_psite_reads to sweep")
    p.add_argument("--paper_frame",  type=float, default=0.55,
                   help="Threshold used in paper (shown as reference line)")
    p.add_argument("--paper_psites", type=int,   default=8,
                   help="Min psites used in paper (shown as reference line)")
    return p.parse_args()

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    df = pd.read_csv(args.input)
    print(f"Loaded {len(df):,} ORFs from {args.input}")

    required = {"total_psites", "frame0_fraction"}
    missing  = required - set(df.columns)
    if missing:
        print(f"ERROR: Missing columns: {missing}")
        print(f"Available: {list(df.columns)}")
        return

    df["total_psites"]    = pd.to_numeric(df["total_psites"],    errors="coerce").fillna(0)
    df["frame0_fraction"] = pd.to_numeric(df["frame0_fraction"], errors="coerce").fillna(0)

    total_orfs = len(df)

    # ── Sweep grid ────────────────────────────────────────────────────────────
    frame_thresholds = np.round(
        np.arange(args.frame_range[0], args.frame_range[1] + 0.01, 0.025), 3
    )
    psite_thresholds = list(range(args.psite_range[0], args.psite_range[1] + 1))

    rows = []
    for ft in frame_thresholds:
        for pt in psite_thresholds:
            passing = df[
                (df["frame0_fraction"] >= ft) &
                (df["total_psites"]    >= pt)
            ]
            n = len(passing)
            rows.append({
                "frame0_threshold": ft,
                "min_psites":       pt,
                "n_passing":        n,
                "pct_passing":      round(100 * n / total_orfs, 2),
            })

    sweep = pd.DataFrame(rows)
    csv_path = os.path.join(args.outdir, f"sensitivity_{args.sample}.csv")
    sweep.to_csv(csv_path, index=False)
    print(f"Sweep table saved: {csv_path}")

    # ── 4-panel figure ────────────────────────────────────────────────────────
    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    fig.suptitle(
        f"Frame-0 threshold sensitivity analysis\n{args.sample}  (N = {total_orfs:,} ORFs)",
        fontsize=13, y=1.01
    )

    # Panel 1: ORF yield vs frame0 threshold (at paper's min_psites)
    ax = axes[0, 0]
    sub = sweep[sweep["min_psites"] == args.paper_psites]
    ax.plot(sub["frame0_threshold"], sub["n_passing"],
            color="steelblue", linewidth=2)
    ax.axvline(args.paper_frame, color="red", linestyle="--",
               label=f"Paper threshold ({args.paper_frame})")
    n_at_paper = sub[sub["frame0_threshold"] == min(
        sub["frame0_threshold"], key=lambda x: abs(x - args.paper_frame)
    )]["n_passing"].values
    if len(n_at_paper):
        ax.axhline(n_at_paper[0], color="red", linestyle=":", alpha=0.5)
    ax.set_xlabel("frame0_fraction threshold")
    ax.set_ylabel("ORFs passing filter")
    ax.set_title(f"ORF yield vs frame0 threshold\n(min_psites = {args.paper_psites})")
    ax.legend(fontsize=9)
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{int(x):,}"))

    # Panel 2: % retained vs frame0 threshold — multiple psite lines
    ax = axes[0, 1]
    psite_levels = [2, 4, 6, 8, 10, 15, 20]
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(psite_levels)))
    for pt, col in zip(psite_levels, colors):
        sub = sweep[sweep["min_psites"] == pt]
        ax.plot(sub["frame0_threshold"], sub["pct_passing"],
                color=col, linewidth=1.5, label=f"min_psites={pt}")
    ax.axvline(args.paper_frame, color="red", linestyle="--",
               linewidth=1.5, label=f"Paper ({args.paper_frame})")
    ax.set_xlabel("frame0_fraction threshold")
    ax.set_ylabel("% ORFs retained")
    ax.set_title("% ORFs retained across thresholds")
    ax.legend(fontsize=7, ncol=2)

    # Panel 3: Heatmap of ORF yield across both thresholds
    ax = axes[1, 0]
    pivot = sweep.pivot(
        index="min_psites", columns="frame0_threshold", values="n_passing"
    )
    im = ax.imshow(
        pivot.values, aspect="auto", cmap="YlOrRd_r",
        origin="lower",
        extent=[frame_thresholds.min(), frame_thresholds.max(),
                psite_thresholds[0], psite_thresholds[-1]]
    )
    ax.axvline(args.paper_frame,  color="cyan",  linestyle="--", linewidth=1.5,
               label=f"paper frame ({args.paper_frame})")
    ax.axhline(args.paper_psites, color="cyan",  linestyle=":",  linewidth=1.5,
               label=f"paper psites ({args.paper_psites})")
    ax.set_xlabel("frame0_fraction threshold")
    ax.set_ylabel("min_psite_reads")
    ax.set_title("ORF yield heatmap")
    ax.legend(fontsize=8)
    plt.colorbar(im, ax=ax, label="ORFs passing")

    # Panel 4: Frame0_fraction distribution of all ORFs
    ax = axes[1, 1]
    ax.hist(df["frame0_fraction"], bins=60, color="steelblue",
            edgecolor="white", linewidth=0.3)
    ax.axvline(args.paper_frame, color="red", linestyle="--",
               linewidth=2, label=f"Paper threshold ({args.paper_frame})")
    ax.axvline(0.333, color="gray", linestyle=":",
               linewidth=1.5, label="Random (0.333)")
    n_above = (df["frame0_fraction"] >= args.paper_frame).sum()
    ax.set_xlabel("frame0_fraction")
    ax.set_ylabel("Number of ORFs")
    ax.set_title(f"Distribution of frame0_fraction\n"
                 f"({n_above:,} ORFs above threshold = "
                 f"{100*n_above/total_orfs:.1f}%)")
    ax.legend(fontsize=9)

    plt.tight_layout()
    pdf_path = os.path.join(args.outdir, f"sensitivity_{args.sample}.pdf")
    plt.savefig(pdf_path, bbox_inches="tight")
    plt.close()
    print(f"Figure saved: {pdf_path}")

    # ── Console summary at paper thresholds ──────────────────────────────────
    at_paper = sweep[
        (sweep["frame0_threshold"] == min(
            sweep["frame0_threshold"].unique(),
            key=lambda x: abs(x - args.paper_frame)
        )) &
        (sweep["min_psites"] == args.paper_psites)
    ]

    print(f"\n{'='*52}")
    print(f"  At paper thresholds (frame0 >= {args.paper_frame}, psites >= {args.paper_psites})")
    if len(at_paper):
        row = at_paper.iloc[0]
        print(f"  ORFs passing : {int(row['n_passing']):,} / {total_orfs:,} "
              f"({row['pct_passing']:.1f}%)")

    # Stable plateau check: find where yield drops < 5% per 0.05 step
    sub_paper_psites = sweep[sweep["min_psites"] == args.paper_psites].sort_values(
        "frame0_threshold"
    )
    yields = sub_paper_psites["n_passing"].values
    thresholds = sub_paper_psites["frame0_threshold"].values
    print(f"\n  Yield stability (min_psites = {args.paper_psites}):")
    for i in range(1, len(yields)):
        if yields[i-1] > 0:
            drop = (yields[i-1] - yields[i]) / yields[i-1] * 100
            marker = " <-- chosen" if abs(thresholds[i] - args.paper_frame) < 0.015 else ""
            print(f"    frame0 >= {thresholds[i]:.3f} : {yields[i]:>6,} ORFs  "
                  f"(drop {drop:+.1f}%){marker}")
    print(f"{'='*52}")

if __name__ == "__main__":
    main()
