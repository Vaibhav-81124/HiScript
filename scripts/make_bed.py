#!/usr/bin/env python3
"""
make_bed.py — Convert cleaned sORF CSV to BED format.

Usage:
    python scripts/make_bed.py \
        --input   results/phase1/stage1_cleaned_sorfs.csv \
        --output  results/phase1/sorfs_genomic.bed
"""

import os
import ast
import argparse
import pandas as pd

def parse_args():
    p = argparse.ArgumentParser(description="sORF CSV → BED")
    p.add_argument("--input",  required=True, help="Cleaned sORF CSV")
    p.add_argument("--output", required=True, help="Output BED file")
    return p.parse_args()

def main():
    args = parse_args()
    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    df   = pd.read_csv(args.input)
    rows = []
    for idx, row in df.iterrows():
        orf_id   = f"ORF_{idx}"
        strand   = row["strand"]
        segments = ast.literal_eval(row["genomic_segments"])
        for chrom, start, end in segments:
            chrom_str = f"chr{chrom}" if not str(chrom).startswith("chr") else str(chrom)
            rows.append([chrom_str, start, end, orf_id, 0, strand])

    bed_df = pd.DataFrame(rows)
    bed_df.to_csv(args.output, sep="\t", header=False, index=False)
    print(f"BED file written: {args.output}  ({len(df):,} ORFs, {len(rows):,} segments)")

if __name__ == "__main__":
    main()
