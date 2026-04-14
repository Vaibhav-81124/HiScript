#!/usr/bin/env python3
"""
merge_reps.py — Intersect translated ORFs across replicates (keep reproducible ORFs).

Usage:
    python scripts/merge_reps.py \
        --rep1    results/phase4/stage5_HeLa_M_rep1_high_confidence_translated_orfs.csv \
        --rep2    results/phase4/stage5_HeLa_M_rep2_high_confidence_translated_orfs.csv \
        --cell    HeLa_M \
        --outdir  results/phase4
"""

import os
import argparse
import pandas as pd

def parse_args():
    p = argparse.ArgumentParser(description="Merge replicate translated ORF tables")
    p.add_argument("--rep1",   required=True, help="Stage 5 CSV rep1")
    p.add_argument("--rep2",   required=True, help="Stage 5 CSV rep2")
    p.add_argument("--cell",   required=True, help="Cell type / condition label")
    p.add_argument("--outdir", required=True, help="Output directory")
    return p.parse_args()

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    rep1 = pd.read_csv(args.rep1)
    rep2 = pd.read_csv(args.rep2)

    print(f"  Rep1 periodic ORFs: {len(rep1):,}")
    print(f"  Rep2 periodic ORFs: {len(rep2):,}")

    common = pd.merge(rep1, rep2, on="genomic_orf_id",
                      suffixes=("_rep1", "_rep2"))

    print(f"  Reproducible ORFs:  {len(common):,}")

    out = os.path.join(args.outdir, f"ribo_{args.cell}_common_translated_orfs.csv")
    common.to_csv(out, index=False)
    print(f"Saved → {out}")

if __name__ == "__main__":
    main()
