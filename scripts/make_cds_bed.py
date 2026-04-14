#!/usr/bin/env python3
"""
make_cds_bed.py — Extract CDS regions from GTF into a BED file.
Run once per genome build. Output is used by stage4_translation.py.

Usage:
    python scripts/make_cds_bed.py \
        --gtf    data/raw/Homo_sapiens.GRCh38.115.gtf \
        --output results/phase1/cds.bed
"""

import os
import re
import argparse

def parse_args():
    p = argparse.ArgumentParser(description="GTF → CDS BED")
    p.add_argument("--gtf",    required=True, help="Genome GTF file")
    p.add_argument("--output", required=True, help="Output CDS BED path")
    return p.parse_args()

def main():
    args = parse_args()
    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    print(f"Extracting CDS from: {args.gtf}")
    count = 0
    with open(args.gtf) as fin, open(args.output, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9 or fields[2] != "CDS":
                continue
            chrom  = fields[0]
            start  = int(fields[3]) - 1      # GTF is 1-based; BED is 0-based
            end    = int(fields[4])
            strand = fields[6]
            m = re.search(r'gene_id "([^"]+)"', fields[8])
            name = m.group(1) if m else "."
            fout.write(f"{chrom}\t{start}\t{end}\t{name}\t0\t{strand}\n")
            count += 1

    print(f"  Written {count:,} CDS intervals → {args.output}")

if __name__ == "__main__":
    main()
