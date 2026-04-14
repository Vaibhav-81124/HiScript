#!/usr/bin/env python3
"""
TE_filter.py — Filter TE table to translated ORFs only; recalculate TE without pseudocount.

Usage:
    python scripts/TE_filter.py \
        --te_table   results/phase5/HeLa_M_translation_efficiency.csv \
        --translated results/phase4/ribo_HeLa_M_common_translated_orfs.csv \
        --cell       HeLa_M \
        --outdir     results/phase5
"""

import os
import argparse
import numpy as np
import pandas as pd

def parse_args():
    p = argparse.ArgumentParser(description="Filter TE to translated ORFs")
    p.add_argument("--te_table",   required=True, help="TE CSV from TE_compute.py")
    p.add_argument("--translated", required=True, help="Common translated ORFs CSV (merge_reps output)")
    p.add_argument("--cell",       required=True, help="Cell type / condition label")
    p.add_argument("--outdir",     required=True, help="Output directory")
    return p.parse_args()

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    df         = pd.read_csv(args.te_table)
    translated = pd.read_csv(args.translated)

    # Keep only reproducibly translated ORFs
    df = df[df["orf_id"].isin(translated["orf_id"])]
    print(f"After translation filtering: {len(df):,}")

    # Recalculate TE without pseudocount (NaN where RNA = 0)
    df["TE_rep1"] = np.where(df["rna_rep1"] > 0,
                             df["ribo_rep1"] / df["rna_rep1"], np.nan)
    df["TE_rep2"] = np.where(df["rna_rep2"] > 0,
                             df["ribo_rep2"] / df["rna_rep2"], np.nan)
    df["TE_mean"]       = df[["TE_rep1","TE_rep2"]].mean(axis=1, skipna=True)
    df["TE_valid_reps"] = df[["TE_rep1","TE_rep2"]].notna().sum(axis=1)

    # Sort by ribosome evidence
    df = df.sort_values(["ribo_rep1","ribo_rep2"], ascending=False)

    out = os.path.join(args.outdir, f"{args.cell}_translated_orfs_filtered_withTE.csv")
    df.to_csv(out, index=False)
    print(f"Final ORFs retained: {len(df):,}")
    print(f"Saved → {out}")

if __name__ == "__main__":
    main()
