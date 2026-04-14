#!/usr/bin/env python3
"""
stage5_periodicity.py — Assess triplet periodicity of Ribo-seq P-sites over candidate ORFs.

Usage:
    python scripts/stage5_periodicity.py \
        --orfs       results/phase4/stage4_HeLa_M_rep1_high_confidence_novel_orfs.csv \
        --psites     results/phase3/HeLa_M_RIBO_rep1_psites.bed \
        --sample     HeLa_M_rep1 \
        --min_reads  8 \
        --frame_thr  0.55 \
        --outdir     results/phase4
"""

import os
import ast
import argparse
import pandas as pd

def parse_args():
    p = argparse.ArgumentParser(description="Stage 5 triplet periodicity filter")
    p.add_argument("--orfs",      required=True, help="Stage 4 high-confidence ORF CSV")
    p.add_argument("--psites",    required=True, help="P-site BED file")
    p.add_argument("--sample",    required=True, help="Sample name")
    p.add_argument("--min_reads", type=int,   default=8,    help="Min P-site reads")
    p.add_argument("--frame_thr", type=float, default=0.55, help="Frame-0 fraction threshold")
    p.add_argument("--outdir",    required=True, help="Output directory")
    return p.parse_args()

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    df = pd.read_csv(args.orfs)
    print(f"Loaded {len(df):,} candidate ORFs")

    psites = pd.read_csv(args.psites, sep="\t", header=None,
                         names=["chr","start","end","name","score","strand"])
    psites["chr"] = psites["chr"].astype(str).str.replace("chr", "", regex=False)

    results = []
    for _, row in df.iterrows():
        segments     = ast.literal_eval(row["genomic_segments"])
        strand       = row["strand"]
        total_reads  = 0
        frame_counts = [0, 0, 0]

        for chrom, start, end in segments:
            seg = psites[
                (psites["chr"]    == str(chrom)) &
                (psites["strand"] == strand) &
                (psites["start"]  >= start) &
                (psites["start"]  <  end)
            ]
            for _, p in seg.iterrows():
                frame = ((p["start"] - start) % 3 if strand == "+"
                         else (end - p["start"] - 1) % 3)
                frame_counts[frame] += 1
                total_reads += 1

        f0_frac = frame_counts[0] / total_reads if total_reads > 0 else 0.0
        results.append([row["orf_id"], total_reads,
                        frame_counts[0], frame_counts[1], frame_counts[2], f0_frac])

    peri_df = pd.DataFrame(results, columns=[
        "orf_id","total_psites","frame0","frame1","frame2","frame0_fraction"])

    merged = df.merge(peri_df, on="orf_id", how="left")
    hc     = merged[
        (merged["total_psites"]    >= args.min_reads) &
        (merged["frame0_fraction"] >= args.frame_thr)
    ]

    print(f"\n{'='*48}")
    print(f"  Total ORFs analyzed:   {len(merged):,}")
    print(f"  Passing periodicity:   {len(hc):,}")
    print(f"  Max frame0_fraction:   {merged['frame0_fraction'].max():.3f}")
    print(f"{'='*48}\n")

    def out(suffix):
        return os.path.join(args.outdir, f"stage5_{args.sample}_{suffix}.csv")

    merged.to_csv(out("with_periodicity"), index=False)
    hc.to_csv(out("high_confidence_translated_orfs"), index=False)
    print(f"Output files written to {args.outdir}/")

if __name__ == "__main__":
    main()
