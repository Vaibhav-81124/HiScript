#!/usr/bin/env python3
"""
sorf_cleaning.py — Phase 1b: Deduplicate and filter raw sORF table.

Usage:
    python scripts/sorf_cleaning.py \
        --input   results/phase1/stage1_novel_sorfs.csv \
        --min_aa  15 \
        --output  results/phase1/stage1_cleaned_sorfs.csv
"""

import os
import argparse
import pandas as pd

def parse_args():
    p = argparse.ArgumentParser(description="Deduplicate and filter sORF table")
    p.add_argument("--input",   required=True, help="Raw sORF CSV from discovery")
    p.add_argument("--min_aa",  type=int, default=15, help="Min aa length after cleaning")
    p.add_argument("--output",  required=True, help="Output cleaned CSV")
    return p.parse_args()

def main():
    args = parse_args()
    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    print(f"Loading: {args.input}")
    df = pd.read_csv(args.input)
    print(f"  Initial ORFs: {len(df):,}")

    # Step 1 — Remove exact duplicate peptides
    df = df.drop_duplicates(subset=["aa_sequence"])
    print(f"  After dedup peptides:    {len(df):,}")

    # Step 2 — Remove duplicate genomic loci
    df = df.drop_duplicates(subset=["chromosome", "strand",
                                    "transcript_start", "transcript_end"])
    print(f"  After dedup loci:        {len(df):,}")

    # Step 3 — Remove nested ORFs (keep longest per chromosome+strand region)
    df["region_key"] = (df["chromosome"].astype(str) + "_" +
                        df["strand"].astype(str))
    kept_rows = []
    for _, group in df.groupby("region_key"):
        group = group.sort_values("aa_length", ascending=False)
        kept  = []
        for _, row in group.iterrows():
            nested = any(
                row["transcript_start"] >= k["transcript_start"] and
                row["transcript_end"]   <= k["transcript_end"]
                for k in kept
            )
            if not nested:
                kept.append(row)
        kept_rows.extend(kept)
    df = pd.DataFrame(kept_rows).drop(columns=["region_key"])
    print(f"  After removing nested:   {len(df):,}")

    # Step 4 — Apply minimum length threshold
    df = df[df["aa_length"] >= args.min_aa]
    print(f"  After min_aa={args.min_aa}:        {len(df):,}")

    df.to_csv(args.output, index=False)
    print(f"\nSaved cleaned sORFs → {args.output}")

if __name__ == "__main__":
    main()
