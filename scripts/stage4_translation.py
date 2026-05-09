#!/usr/bin/env python3
"""
stage4_translation.py — Merge Ribo counts, remove CDS overlap, filter by read threshold.

Now uses stable genomic_orf_id across the entire pipeline.
"""

import os
import ast
import hashlib
import subprocess
import argparse
import pandas as pd

def parse_args():
    p = argparse.ArgumentParser(description="Stage 4 translation filter")
    p.add_argument("--sorfs",     required=True, help="Cleaned sORF CSV")
    p.add_argument("--ribo",      required=True, help="Ribo count txt (genomic_orf_id  count)")
    p.add_argument("--cds_bed",   required=True, help="CDS BED file")
    p.add_argument("--sorf_bed",  required=True, help="sORF genomic BED file")
    p.add_argument("--sample",    required=True, help="Sample name")
    p.add_argument("--min_reads", type=int, default=10)
    p.add_argument("--micro_max", type=int, default=60)
    p.add_argument("--outdir",    required=True)
    return p.parse_args()

# ─────────────────────────────────────────────────────────────

def make_genomic_id(row):
    segments = ast.literal_eval(row["genomic_segments"]) if isinstance(row["genomic_segments"], str) else row["genomic_segments"]
    norm_segments = []
    for c, s, e in segments:
        # Strip chr prefix - consistent with Ensembl no-chr convention
        chrom = str(c).replace("chr", "")
        norm_segments.append((chrom, s, e))

    seg_str = ";".join([f"{c}:{s}-{e}" for c, s, e in norm_segments])
    return f"{seg_str}|{row['strand']}"

def hash_peptide(seq):
    return hashlib.md5(seq.encode()).hexdigest()[:10]

# ─────────────────────────────────────────────────────────────

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # ── Load sORFs ──────────────────────────────────────────
    df = pd.read_csv(args.sorfs)
    print(f"Loaded {len(df):,} sORFs")

    # ── Create stable IDs ───────────────────────────────────
    df["genomic_orf_id"] = df.apply(make_genomic_id, axis=1)
    df["peptide_orf_id"] = df["aa_sequence"].apply(hash_peptide)

    # ── Merge Ribo counts ───────────────────────────────────
    counts = pd.read_csv(
        args.ribo,
        sep=r"\s+",
        header=None,
        names=["genomic_orf_id", "ribo_count"]
    )

    df = df.merge(counts, on="genomic_orf_id", how="left")
    df["ribo_count"] = pd.to_numeric(df["ribo_count"], errors="coerce").fillna(0).astype(int)

    print("Ribo counts merged.")

    # 🔴 Sanity check (prevents silent failure)
    nonzero = (df["ribo_count"] > 0).sum()
    print(f"Non-zero ribo counts: {nonzero:,}")
    if nonzero == 0:
        raise ValueError("All ribo counts are zero → ID mismatch between BED and count file.")

    # ── CDS overlap removal ─────────────────────────────────
    print("Removing CDS-overlapping ORFs...")

    tmp = os.path.join(args.outdir, "temp_cds_overlap.bed")

    subprocess.run(
        f"bedtools intersect -a {args.sorf_bed} -b {args.cds_bed} -u > {tmp}",
        shell=True,
        check=True
    )

    remove_cds = set(pd.read_csv(tmp, sep="\t", header=None)[3])

    df_novel = df[~df["genomic_orf_id"].isin(remove_cds)]
    os.remove(tmp)

    print(f"  After CDS removal: {len(df_novel):,}")

    # ── Translation filter ──────────────────────────────────
    df_translated = df_novel[df_novel["ribo_count"] >= args.min_reads]
    df_micro      = df_translated[df_translated["aa_length"] <= args.micro_max]

    # ── Summary ─────────────────────────────────────────────
    print(f"\n{'='*50}")
    print(f"  Total ORFs:                        {len(df):,}")
    print(f"  After CDS removal:                 {len(df_novel):,}")
    print(f"  High-confidence (≥{args.min_reads}): {len(df_translated):,}")
    print(f"  Microproteins (≤{args.micro_max} aa): {len(df_micro):,}")
    print(f"{'='*50}\n")

    # ── Save ────────────────────────────────────────────────
    def out(suffix):
        return os.path.join(args.outdir, f"stage4_{args.sample}_{suffix}.csv")

    df.to_csv(out("all_with_ribo"), index=False)
    df_novel.to_csv(out("novel_orfs"), index=False)
    df_translated.to_csv(out("high_confidence_novel_orfs"), index=False)
    df_micro.to_csv(out("high_confidence_microproteins"), index=False)

    print(f"Output files written to {args.outdir}/")

# ─────────────────────────────────────────────────────────────

if __name__ == "__main__":
    main()