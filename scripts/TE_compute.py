#!/usr/bin/env python3
"""
TE_compute.py — Compute Translation Efficiency (TE) from Ribo and RNA count tables.

Usage:
    python scripts/TE_compute.py \
        --ribo1   results/phase4/sorf_total_ribo_counts_HeLa_M_rep1.txt \
        --ribo2   results/phase4/sorf_total_ribo_counts_HeLa_M_rep2.txt \
        --rna1    results/phase2/rna_count/sorf_rna_counts_HeLa_M_rep1_clean.txt \
        --rna2    results/phase2/rna_count/sorf_rna_counts_HeLa_M_rep2_clean.txt \
        --cell    HeLa_M \
        --pseudo  1 \
        --outdir  results/phase5
"""

import os
import argparse
import pandas as pd

def parse_args():
    p = argparse.ArgumentParser(description="Compute Translation Efficiency")
    p.add_argument("--ribo1",   required=True)
    p.add_argument("--ribo2",   required=True)
    p.add_argument("--rna1",    required=True)
    p.add_argument("--rna2",    required=True)
    p.add_argument("--cell",    required=True, help="Cell type / condition label")
    p.add_argument("--pseudo",  type=int, default=1, help="Pseudocount added to RNA denominator")
    p.add_argument("--outdir",  required=True)
    return p.parse_args()

def load_counts(path, col_name):
    return pd.read_csv(path, sep=r"\s+", header=None,
                       names=["orf_id", col_name])

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    ribo1 = load_counts(args.ribo1, "ribo_rep1")
    ribo2 = load_counts(args.ribo2, "ribo_rep2")
    rna1  = load_counts(args.rna1,  "rna_rep1")
    rna2  = load_counts(args.rna2,  "rna_rep2")

    df = ribo1.merge(ribo2, on="orf_id", how="outer") \
              .merge(rna1,  on="orf_id", how="outer") \
              .merge(rna2,  on="orf_id", how="outer")
    df = df.fillna(0)

    for col in ["ribo_rep1","ribo_rep2","rna_rep1","rna_rep2"]:
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0)

    df["TE_rep1"] = df["ribo_rep1"] / (df["rna_rep1"] + args.pseudo)
    df["TE_rep2"] = df["ribo_rep2"] / (df["rna_rep2"] + args.pseudo)
    df["TE_mean"] = (df["TE_rep1"] + df["TE_rep2"]) / 2

    out = os.path.join(args.outdir, f"{args.cell}_translation_efficiency.csv")
    df.to_csv(out, index=False)
    print(f"TE table written → {out}  ({len(df):,} ORFs)")

if __name__ == "__main__":
    main()
