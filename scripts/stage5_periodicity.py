#!/usr/bin/env python3
"""
stage5_periodicity.py — Assess triplet periodicity of Ribo-seq P-sites over candidate ORFs.

Now uses genomic_orf_id (stable ID) instead of index-based ORF IDs.
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
    p.add_argument("--min_reads", type=int,   default=8)
    p.add_argument("--frame_thr", type=float, default=0.55)
    p.add_argument("--outdir",    required=True)
    return p.parse_args()

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # ── Load ORFs ──────────────────────────────────────────
    df = pd.read_csv(args.orfs)
    print(f"Loaded {len(df):,} candidate ORFs")

    if "genomic_orf_id" not in df.columns:
        raise ValueError("Missing genomic_orf_id in input ORF file. Fix Stage 4 output.")

    # ── Load P-sites ───────────────────────────────────────
    psites = pd.read_csv(
        args.psites,
        sep="\t",
        header=None,
        names=["chr","start","end","name","score","strand"]
    )

    # Normalize chromosome format
    psites["chr"] = psites["chr"].astype(str).str.replace("chr", "", regex=False)

    # ── Compute periodicity ────────────────────────────────
    results = []

    for _, row in df.iterrows():
        segments = ast.literal_eval(row["genomic_segments"]) \
            if isinstance(row["genomic_segments"], str) else row["genomic_segments"]

        strand = row["strand"]

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
                if strand == "+":
                    frame = (p["start"] - start) % 3
                else:
                    frame = (end - p["start"] - 1) % 3

                frame_counts[frame] += 1
                total_reads += 1

        f0_frac = frame_counts[0] / total_reads if total_reads > 0 else 0.0

        results.append([
            row["genomic_orf_id"],
            total_reads,
            frame_counts[0],
            frame_counts[1],
            frame_counts[2],
            f0_frac
        ])

    # ── Build dataframe ────────────────────────────────────
    peri_df = pd.DataFrame(results, columns=[
        "genomic_orf_id",
        "total_psites",
        "frame0",
        "frame1",
        "frame2",
        "frame0_fraction"
    ])

    # ── Merge back ─────────────────────────────────────────
    merged = df.merge(peri_df, on="genomic_orf_id", how="left")

    # ── Apply filters ──────────────────────────────────────
    hc = merged[
        (merged["total_psites"]    >= args.min_reads) &
        (merged["frame0_fraction"] >= args.frame_thr)
    ]

    # ── Summary ────────────────────────────────────────────
    print(f"\n{'='*50}")
    print(f"  Total ORFs analyzed:   {len(merged):,}")
    print(f"  Passing periodicity:   {len(hc):,}")
    print(f"  Max frame0_fraction:   {merged['frame0_fraction'].max():.3f}")
    print(f"{'='*50}\n")

    # ── Save ───────────────────────────────────────────────
    def out(suffix):
        return os.path.join(args.outdir, f"stage5_{args.sample}_{suffix}.csv")

    merged.to_csv(out("with_periodicity"), index=False)
    hc.to_csv(out("high_confidence_translated_orfs"), index=False)

    print(f"Output files written to {args.outdir}/")

if __name__ == "__main__":
    main()