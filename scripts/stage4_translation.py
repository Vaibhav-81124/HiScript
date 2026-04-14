#!/usr/bin/env python3
"""
stage4_translation.py — Merge Ribo counts, remove CDS overlap, filter by read threshold.

Usage:
    python scripts/stage4_translation.py \
        --sorfs      results/phase1/stage1_cleaned_sorfs.csv \
        --ribo       results/phase4/sorf_total_ribo_counts_rep1.txt \
        --cds_bed    results/phase1/cds.bed \
        --sorf_bed   results/phase1/sorfs_genomic.bed \
        --sample     HeLa_M_rep1 \
        --min_reads  10 \
        --micro_max  60 \
        --outdir     results/phase4
"""

import os
import hashlib
import subprocess
import argparse
import pandas as pd

def parse_args():
    p = argparse.ArgumentParser(description="Stage 4 translation filter")
    p.add_argument("--sorfs",     required=True, help="Cleaned sORF CSV")
    p.add_argument("--ribo",      required=True, help="Ribo count txt (orf_id  count)")
    p.add_argument("--cds_bed",   required=True, help="CDS BED file")
    p.add_argument("--sorf_bed",  required=True, help="sORF genomic BED file")
    p.add_argument("--sample",    required=True, help="Sample name (used in output filenames)")
    p.add_argument("--min_reads", type=int, default=10, help="Min Ribo reads to call translated")
    p.add_argument("--micro_max", type=int, default=60, help="Max aa for microprotein class")
    p.add_argument("--outdir",    required=True, help="Output directory")
    return p.parse_args()

def make_genomic_id(row):
    return f"{row['genomic_segments']}_{row['strand']}"

def hash_peptide(seq):
    return hashlib.md5(seq.encode()).hexdigest()[:10]

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # ── Load sORFs ──────────────────────────────────────────────────────────
    df = pd.read_csv(args.sorfs)
    print(f"Loaded {len(df):,} sORFs")
    df["orf_id"] = df.index.map(lambda x: f"ORF_{x}")

    # ── Merge Ribo counts ───────────────────────────────────────────────────
    counts = pd.read_csv(args.ribo, sep=r"\s+", header=None,
                         names=["orf_id", "ribo_count"])
    df = df.merge(counts, on="orf_id", how="left")
    df["ribo_count"] = pd.to_numeric(df["ribo_count"], errors="coerce").fillna(0).astype(int)
    print("Ribo counts merged.")

    # ── Stable IDs ──────────────────────────────────────────────────────────
    df["genomic_orf_id"] = df.apply(make_genomic_id, axis=1)
    df["peptide_orf_id"] = df["aa_sequence"].apply(hash_peptide)

    # ── CDS overlap removal ─────────────────────────────────────────────────
    print("Removing CDS-overlapping ORFs...")
    tmp = os.path.join(args.outdir, "temp_cds_overlap.bed")
    subprocess.run(
        f"bedtools intersect -a {args.sorf_bed} -b {args.cds_bed} -u > {tmp}",
        shell=True, check=True
    )
    remove_cds = set(pd.read_csv(tmp, sep="\t", header=None)[3])
    df_novel   = df[~df["orf_id"].isin(remove_cds)]
    os.remove(tmp)
    print(f"  After CDS removal: {len(df_novel):,}")

    # ── Translation filter ──────────────────────────────────────────────────
    df_translated = df_novel[df_novel["ribo_count"] >= args.min_reads]
    df_micro      = df_translated[df_translated["aa_length"] <= args.micro_max]

    # ── Summary ─────────────────────────────────────────────────────────────
    print(f"\n{'='*48}")
    print(f"  Total ORFs:                    {len(df):,}")
    print(f"  After CDS removal:             {len(df_novel):,}")
    print(f"  High-confidence (≥{args.min_reads} reads): {len(df_translated):,}")
    print(f"  Microproteins (≤{args.micro_max} aa):     {len(df_micro):,}")
    print(f"{'='*48}\n")

    # ── Save ─────────────────────────────────────────────────────────────────
    def out(suffix):
        return os.path.join(args.outdir, f"stage4_{args.sample}_{suffix}.csv")

    df.to_csv(out("all_with_ribo"), index=False)
    df_novel.to_csv(out("novel_orfs"), index=False)
    df_translated.to_csv(out("high_confidence_novel_orfs"), index=False)
    df_micro.to_csv(out("high_confidence_microproteins"), index=False)
    print(f"Output files written to {args.outdir}/")

if __name__ == "__main__":
    main()
