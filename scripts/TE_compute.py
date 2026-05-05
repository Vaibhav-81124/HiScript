#!/usr/bin/env python3
"""
TE_compute.py — Compute Translation Efficiency (TE) from Ribo and RNA count tables.

Now uses genomic_orf_id (stable ID) instead of index-based ORF IDs.
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
    df = pd.read_csv(path, sep=r"\s+", header=None)

    if df.shape[1] < 2:
        raise ValueError(f"{path} must have at least 2 columns")

    df = df.iloc[:, :2]
    df.columns = ["genomic_orf_id", col_name]

    return df

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # ── Load counts ───────────────────────────────────────
    ribo1 = load_counts(args.ribo1, "ribo_rep1")
    ribo2 = load_counts(args.ribo2, "ribo_rep2")
    rna1  = load_counts(args.rna1,  "rna_rep1")
    rna2  = load_counts(args.rna2,  "rna_rep2")

    # ── Merge ─────────────────────────────────────────────
    df = ribo1.merge(ribo2, on="genomic_orf_id", how="outer") \
              .merge(rna1,  on="genomic_orf_id", how="outer") \
              .merge(rna2,  on="genomic_orf_id", how="outer")

    df = df.fillna(0)

    # ── Convert to numeric ────────────────────────────────
    for col in ["ribo_rep1","ribo_rep2","rna_rep1","rna_rep2"]:
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0)

    # ── Sanity check ──────────────────────────────────────
    nonzero_ribo = (df["ribo_rep1"] > 0).sum() + (df["ribo_rep2"] > 0).sum()
    print(f"Non-zero ribo entries: {nonzero_ribo:,}")

    if nonzero_ribo == 0:
        raise ValueError(
            "All Ribo counts are zero after merge.\n"
            "Likely ID mismatch between ribo/RNA files."
        )

    # ── Compute TE ────────────────────────────────────────
    df["TE_rep1"] = df["ribo_rep1"] / (df["rna_rep1"] + args.pseudo)
    df["TE_rep2"] = df["ribo_rep2"] / (df["rna_rep2"] + args.pseudo)
    df["TE_mean"] = (df["TE_rep1"] + df["TE_rep2"]) / 2

    # ── Save ──────────────────────────────────────────────
    out = os.path.join(args.outdir, f"{args.cell}_translation_efficiency.csv")
    df.to_csv(out, index=False)

    print(f"TE table written → {out}  ({len(df):,} ORFs)")

if __name__ == "__main__":
    main()