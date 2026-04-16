#!/usr/bin/env python3
"""
psite_calibration.py — Empirically determine optimal P-site offsets from your own data.

HOW IT WORKS:
  Annotated CDS regions have a known reading frame. A correctly assigned P-site
  should land in frame 0 with respect to the CDS start codon. For each RPF
  read length (26–34 nt), we test every possible offset (10–16) and pick the
  one that maximises the fraction of reads landing in frame 0 over all annotated
  CDSs. These empirical offsets are then written to a JSON file that
  psite_assignment.py (which replaces the awk step in run_phase3.sh) reads at
  runtime — so every sample self-calibrates from its own BAM.

PIPELINE POSITION:
  Run AFTER STAR alignment + RPF BAM filtering (Phase 3 step 4),
  BEFORE P-site assignment (Phase 3 step 5).

  Phase 3 modified order:
    trim → rRNA remove → STAR → RPF filter → [psite_calibration] → psite_assignment

Usage:
    python scripts/psite_calibration.py \
        --bam        results/phase3/03_aligned/HeLa_M_RIBO_rep1_RPFs.bam \
        --gtf        data/raw/Homo_sapiens.GRCh38.115.gtf \
        --sample     HeLa_M_RIBO_rep1 \
        --rpf_min    26 \
        --rpf_max    34 \
        --outdir     results/phase3/04_psites
"""

import os
import re
import json
import argparse
import subprocess
import tempfile
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description="Empirical P-site offset calibration")
    p.add_argument("--bam",      required=True,  help="RPF-filtered BAM file")
    p.add_argument("--gtf",      required=True,  help="Genome GTF (for CDS annotation)")
    p.add_argument("--sample",   required=True,  help="Sample name")
    p.add_argument("--rpf_min",  type=int, default=26)
    p.add_argument("--rpf_max",  type=int, default=34)
    p.add_argument("--min_cds_len",  type=int, default=100,
                   help="Min CDS length (nt) to use for calibration")
    p.add_argument("--n_cds",    type=int, default=5000,
                   help="Number of CDS regions to use (top expressed)")
    p.add_argument("--outdir",   required=True,  help="Output directory")
    return p.parse_args()

# ── PARSE CDS FROM GTF ────────────────────────────────────────────────────────

def parse_cds_regions(gtf_file: str, min_len: int):
    """
    Extract CDS start codon positions from GTF.
    Returns list of (chrom, cds_start_genome, strand, gene_id).
    We use just the first 90 nt of each CDS (start codon region) for calibration
    because that's where frame signal is strongest and cleanest.
    """
    print(f"Parsing CDS regions from GTF...")
    cds_dict = defaultdict(list)  # gene_id → list of (chrom, start, end, strand)

    with open(gtf_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9 or fields[2] != "CDS":
                continue
            chrom  = fields[0]
            start  = int(fields[3]) - 1   # 0-based
            end    = int(fields[4])
            strand = fields[6]
            m = re.search(r'gene_id "([^"]+)"', fields[8])
            if not m:
                continue
            gene_id = m.group(1)
            cds_dict[gene_id].append((chrom, start, end, strand))

    # Flatten to one entry per gene: full CDS span
    regions = []
    for gene_id, intervals in cds_dict.items():
        chrom  = intervals[0][0]
        strand = intervals[0][3]
        g_min  = min(s for _, s, _, _ in intervals)
        g_max  = max(e for _, _, e, _ in intervals)
        length = g_max - g_min
        if length >= min_len:
            regions.append({
                "gene_id": gene_id,
                "chrom":   chrom,
                "start":   g_min,
                "end":     g_max,
                "strand":  strand,
                "length":  length,
            })

    df = pd.DataFrame(regions)
    print(f"  {len(df):,} CDS regions found (>= {min_len} nt)")
    return df

# ── EXTRACT READS OVER CDS VIA SAMTOOLS ──────────────────────────────────────

def get_reads_over_cds(bam: str, cds_df: pd.DataFrame,
                        rpf_min: int, rpf_max: int, n_cds: int):
    """
    For each CDS region, extract aligned reads using samtools view.
    Returns a list of dicts: {chrom, pos, strand, readlen, cds_start, cds_strand}
    """
    print(f"Extracting reads over {min(n_cds, len(cds_df)):,} CDS regions...")

    # Sample top n_cds regions (longest = best expressed proxy)
    sample = cds_df.nlargest(n_cds, "length")

    reads = []
    for _, row in sample.iterrows():
        region = f"{row['chrom']}:{row['start']+1}-{row['end']}"
        try:
            result = subprocess.run(
                ["samtools", "view", bam, region],
                capture_output=True, text=True, check=True
            )
        except subprocess.CalledProcessError:
            continue

        for line in result.stdout.strip().split("\n"):
            if not line:
                continue
            fields  = line.split("\t")
            if len(fields) < 10:
                continue
            flag    = int(fields[1])
            pos     = int(fields[3]) - 1    # 0-based
            readlen = len(fields[9])
            strand  = "-" if (flag & 16) else "+"

            if rpf_min <= readlen <= rpf_max:
                reads.append({
                    "chrom":      row["chrom"],
                    "pos":        pos,
                    "strand":     strand,
                    "readlen":    readlen,
                    "cds_start":  row["start"],
                    "cds_end":    row["end"],
                    "cds_strand": row["strand"],
                })

    print(f"  {len(reads):,} reads extracted")
    return reads

# ── CALIBRATE OFFSETS ────────────────────────────────────────────────────────

def calibrate_offsets(reads: list, rpf_min: int, rpf_max: int):
    """
    For each read length, test offsets 10–16.
    Score = fraction of reads landing in frame 0 of the CDS.
    Best offset = argmax(frame0_fraction).
    """
    print("Calibrating offsets...")

    # Group reads by read length
    by_len = defaultdict(list)
    for r in reads:
        by_len[r["readlen"]].append(r)

    offset_candidates = list(range(10, 17))
    results = {}

    for rlen in range(rpf_min, rpf_max + 1):
        rlen_reads = by_len.get(rlen, [])
        if len(rlen_reads) < 50:
            # Not enough reads — use default offset
            results[rlen] = {"best_offset": 12, "frame0_frac": None,
                             "n_reads": len(rlen_reads), "calibrated": False}
            continue

        best_offset    = 12
        best_f0_frac   = 0.0
        offset_scores  = {}

        for offset in offset_candidates:
            frame0_count = 0
            total        = 0

            for r in rlen_reads:
                # Only use reads on same strand as CDS
                if r["strand"] != r["cds_strand"]:
                    continue

                if r["strand"] == "+":
                    psite = r["pos"] + offset
                    frame = (psite - r["cds_start"]) % 3
                else:
                    psite = r["pos"] + rlen - offset - 1
                    frame = (r["cds_end"] - psite - 1) % 3

                if frame == 0:
                    frame0_count += 1
                total += 1

            f0_frac = frame0_count / total if total > 0 else 0.0
            offset_scores[offset] = f0_frac

            if f0_frac > best_f0_frac:
                best_f0_frac = f0_frac
                best_offset  = offset

        results[rlen] = {
            "best_offset":    best_offset,
            "frame0_frac":    round(best_f0_frac, 4),
            "offset_scores":  {str(k): round(v, 4) for k, v in offset_scores.items()},
            "n_reads":        len(rlen_reads),
            "calibrated":     True,
        }

        print(f"  Read length {rlen}: best offset = {best_offset}  "
              f"(frame0 = {best_f0_frac:.3f}, n = {len(rlen_reads):,})")

    return results

# ── PLOT ─────────────────────────────────────────────────────────────────────

def plot_calibration(results: dict, sample: str, outdir: str):
    lengths   = sorted(results.keys())
    offsets   = [results[l]["best_offset"] for l in lengths]
    f0_fracs  = [results[l]["frame0_frac"] if results[l]["frame0_frac"] else 0
                 for l in lengths]
    calibrated = [results[l]["calibrated"] for l in lengths]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 7), sharex=True)

    colors = ["steelblue" if c else "lightgray" for c in calibrated]

    ax1.bar(lengths, offsets, color=colors, edgecolor="white", linewidth=0.5)
    ax1.axhline(12, color="red", linestyle="--", linewidth=1, alpha=0.6, label="Default offset (12)")
    ax1.set_ylabel("Best P-site offset (nt)")
    ax1.set_title(f"{sample} — Empirical P-site offsets per read length")
    ax1.legend(fontsize=9)
    ax1.set_ylim(9, 17)
    ax1.set_yticks(range(10, 17))

    ax2.bar(lengths, f0_fracs, color=colors, edgecolor="white", linewidth=0.5)
    ax2.axhline(0.333, color="red", linestyle="--", linewidth=1, alpha=0.6, label="Random (0.333)")
    ax2.set_ylabel("Frame-0 fraction")
    ax2.set_xlabel("RPF read length (nt)")
    ax2.legend(fontsize=9)
    ax2.set_ylim(0, 1.0)

    # Add legend for calibration status
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor="steelblue", label="Calibrated"),
        Patch(facecolor="lightgray", label="Default (insufficient reads)")
    ]
    ax1.legend(handles=legend_elements + [
        plt.Line2D([0], [0], color="red", linestyle="--", label="Default offset (12)")
    ], fontsize=8)

    plt.tight_layout()
    plot_path = os.path.join(outdir, f"{sample}_psite_calibration.pdf")
    plt.savefig(plot_path, bbox_inches="tight")
    plt.close()
    print(f"  Plot saved: {plot_path}")

# ── MAIN ──────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # Step 1: Parse CDS
    cds_df = parse_cds_regions(args.gtf, args.min_cds_len)

    # Step 2: Extract reads over CDS
    reads = get_reads_over_cds(
        args.bam, cds_df, args.rpf_min, args.rpf_max, args.n_cds
    )

    if len(reads) < 500:
        print(f"WARNING: Only {len(reads)} reads found over CDS regions.")
        print("Calibration may be unreliable. Check BAM file and CDS overlap.")

    # Step 3: Calibrate
    results = calibrate_offsets(reads, args.rpf_min, args.rpf_max)

    # Step 4: Save offsets JSON
    offsets_simple = {str(rlen): v["best_offset"] for rlen, v in results.items()}
    json_path = os.path.join(args.outdir, f"{args.sample}_psite_offsets.json")
    with open(json_path, "w") as f:
        json.dump(offsets_simple, f, indent=2)
    print(f"\nOffset table saved: {json_path}")

    # Step 5: Save full calibration report
    report_rows = []
    for rlen, v in sorted(results.items()):
        report_rows.append({
            "read_length":  rlen,
            "best_offset":  v["best_offset"],
            "frame0_frac":  v["frame0_frac"],
            "n_reads":      v["n_reads"],
            "calibrated":   v["calibrated"],
        })
    report_df = pd.DataFrame(report_rows)
    report_path = os.path.join(args.outdir, f"{args.sample}_calibration_report.csv")
    report_df.to_csv(report_path, index=False)
    print(f"Calibration report saved: {report_path}")

    # Step 6: Plot
    try:
        plot_calibration(results, args.sample, args.outdir)
    except Exception as e:
        print(f"  (Plot skipped: {e})")

    # Step 7: Summary
    calibrated_lens = [r for r in results.values() if r["calibrated"]]
    default_lens    = [r for r in results.values() if not r["calibrated"]]
    print(f"\n{'='*52}")
    print(f"  Calibration Summary for {args.sample}")
    print(f"{'='*52}")
    print(f"  Read lengths calibrated : {len(calibrated_lens)}")
    print(f"  Using default offset    : {len(default_lens)} (insufficient reads)")
    if calibrated_lens:
        mean_f0 = np.mean([r["frame0_frac"] for r in calibrated_lens])
        print(f"  Mean frame-0 fraction   : {mean_f0:.3f}  (random = 0.333)")
        if mean_f0 > 0.50:
            print(f"  ✓ Strong periodicity signal detected")
        elif mean_f0 > 0.40:
            print(f"  ~ Moderate periodicity signal")
        else:
            print(f"  ✗ Weak periodicity — check library quality")
    print(f"{'='*52}")
    print(f"\nNext step: pass {json_path} to psite_assignment.py")

if __name__ == "__main__":
    main()
