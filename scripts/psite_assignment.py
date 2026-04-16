#!/usr/bin/env python3
"""
psite_assignment.py — Assign P-sites using empirically calibrated offsets.

Replaces the hardcoded awk P-site step in run_phase3.sh.
Reads the offset JSON produced by psite_calibration.py — if no JSON is
provided, falls back to the original hardcoded offsets so it's backward
compatible.

Usage:
    python scripts/psite_assignment.py \
        --bam      results/phase3/03_aligned/HeLa_M_RIBO_rep1_RPFs.bam \
        --offsets  results/phase3/04_psites/HeLa_M_RIBO_rep1_psite_offsets.json \
        --sample   HeLa_M_RIBO_rep1 \
        --rpf_min  26 \
        --rpf_max  34 \
        --outdir   results/phase3/04_psites

    # Fallback (no offsets JSON — uses default table):
    python scripts/psite_assignment.py \
        --bam      results/phase3/03_aligned/HeLa_M_RIBO_rep1_RPFs.bam \
        --sample   HeLa_M_RIBO_rep1 \
        --outdir   results/phase3/04_psites
"""

import os
import json
import argparse
import subprocess

# ── Default hardcoded offsets (original pipeline) ────────────────────────────
DEFAULT_OFFSETS = {
    26: 12, 27: 12, 28: 12, 29: 12,
    30: 13, 31: 13,
    32: 12, 33: 12, 34: 12,
}

def parse_args():
    p = argparse.ArgumentParser(description="P-site assignment with calibrated offsets")
    p.add_argument("--bam",     required=True,  help="RPF-filtered BAM")
    p.add_argument("--offsets", default=None,   help="Offsets JSON from psite_calibration.py")
    p.add_argument("--sample",  required=True)
    p.add_argument("--rpf_min", type=int, default=26)
    p.add_argument("--rpf_max", type=int, default=34)
    p.add_argument("--outdir",  required=True)
    return p.parse_args()

def load_offsets(json_path, rpf_min, rpf_max):
    if json_path and os.path.exists(json_path):
        with open(json_path) as f:
            raw = json.load(f)
        offsets = {int(k): int(v) for k, v in raw.items()}
        print(f"  Using calibrated offsets from: {json_path}")
        source = "calibrated"
    else:
        offsets = DEFAULT_OFFSETS.copy()
        print(f"  No offsets JSON found — using default offsets")
        source = "default"

    # Fill any missing lengths with 12
    for l in range(rpf_min, rpf_max + 1):
        if l not in offsets:
            offsets[l] = 12

    # Print the table
    print(f"\n  {'Length':>8}  {'Offset':>8}  {'Source':>12}")
    print(f"  {'─'*8}  {'─'*8}  {'─'*12}")
    for l in sorted(offsets.keys()):
        src = source if l in offsets else "default"
        print(f"  {l:>8}  {offsets[l]:>8}  {src:>12}")
    print()

    return offsets

def assign_psites(bam, offsets, rpf_min, rpf_max, outdir, sample):
    """
    Stream through BAM with samtools view and assign P-sites.
    Writes BED6 output: chrom  psite  psite+1  .  .  strand
    """
    out_path = os.path.join(outdir, f"{sample}_psites.bed")
    print(f"  Assigning P-sites → {out_path}")

    proc = subprocess.Popen(
        ["samtools", "view", bam],
        stdout=subprocess.PIPE,
        text=True
    )

    total   = 0
    written = 0

    with open(out_path, "w") as fout:
        for line in proc.stdout:
            fields  = line.split("\t")
            if len(fields) < 10:
                continue
            flag    = int(fields[1])
            chrom   = fields[2]
            pos     = int(fields[3]) - 1    # 0-based
            readlen = len(fields[9])
            strand  = "-" if (flag & 16) else "+"

            total += 1
            if readlen < rpf_min or readlen > rpf_max:
                continue

            offset = offsets.get(readlen, 12)

            if strand == "+":
                psite = pos + offset
            else:
                psite = pos + readlen - offset - 1

            fout.write(f"{chrom}\t{psite}\t{psite+1}\t.\t.\t{strand}\n")
            written += 1

    proc.wait()
    print(f"  Total reads processed : {total:,}")
    print(f"  P-sites written       : {written:,}")
    return out_path

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    print(f"\nP-site assignment: {args.sample}")
    offsets  = load_offsets(args.offsets, args.rpf_min, args.rpf_max)
    out_path = assign_psites(
        args.bam, offsets, args.rpf_min, args.rpf_max,
        args.outdir, args.sample
    )
    print(f"\n✓ P-site BED written: {out_path}")

if __name__ == "__main__":
    main()
