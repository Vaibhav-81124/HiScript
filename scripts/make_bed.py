#!/usr/bin/env python3
"""
make_bed.py — Convert cleaned sORF CSV to BED format using stable genomic IDs.

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
    p = argparse.ArgumentParser(description="sORF CSV → BED (stable IDs)")
    p.add_argument("--input",  required=True, help="Cleaned sORF CSV")
    p.add_argument("--output", required=True, help="Output BED file")
    return p.parse_args()


def make_genomic_id(segments, strand):
    """
    Create a stable genomic ORF ID.

    Format:
        chr:start-end;chr:start-end|strand

    Example:
        chr1:100-150;chr1:200-250|+
    """
    seg_str = ";".join([f"{c}:{s}-{e}" for c, s, e in segments])
    return f"{seg_str}|{strand}"


def main():
    args = parse_args()
    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    df = pd.read_csv(args.input)
    rows = []

    for _, row in df.iterrows():
        strand = row["strand"]

        # Parse segments safely
        segments = row["genomic_segments"]
        if isinstance(segments, str):
            segments = ast.literal_eval(segments)

        # Normalize chromosome naming - strip chr prefix for consistency
        # STAR with Ensembl GTF outputs chroms without chr prefix (1, 2, X etc.)
        norm_segments = []
        for chrom, start, end in segments:
            chrom_str = str(chrom).replace("chr", "")
            norm_segments.append((chrom_str, int(start), int(end)))

        # Stable genomic ID
        genomic_id = make_genomic_id(norm_segments, strand)

        # Write each segment as BED row
        for chrom_str, start, end in norm_segments:
            rows.append([chrom_str, start, end, genomic_id, 0, strand])

    bed_df = pd.DataFrame(rows)

    bed_df.to_csv(args.output, sep="\t", header=False, index=False)

    print(f"BED file written: {args.output}")
    print(f"  ORFs:     {len(df):,}")
    print(f"  Segments: {len(rows):,}")


if __name__ == "__main__":
    main()