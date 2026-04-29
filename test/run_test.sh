#!/usr/bin/env bash
# =============================================================================
#  run_test.sh — Quick validation run using public HeLa data (SRA)
#
#  Uses two HeLa M-phase replicates (RNA-seq + Ribo-seq) from:
#    Aviner et al. 2017, Nat. Struct. Mol. Biol.
#    SRR3306577, SRR3306578  (RNA-seq)
#    SRR3306585, SRR3306586  (Ribo-seq)
#
#  What this test does:
#    - Downloads ~500MB of FASTQ data from SRA automatically
#    - Skips Phase 1 (uses precomputed cleaned sORF table from Releases)
#    - Runs Phases 2–5 on the HeLa_M samples
#    - Full run takes ~1–2 hours depending on your machine
#
#  Prerequisites:
#    1. conda env create -f environment.yml && conda activate sorf-tool
#    2. Place GRCh38 reference files in data/raw/ (see README)
#    3. Download stage1_cleaned_sorfs.csv from GitHub Releases into results/phase1/
#
#  Usage:
#    bash test/run_test.sh
# =============================================================================

set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT}"

echo "╔══════════════════════════════════════════════════════════╗"
echo "║           sORF Tool — Test Run (HeLa M-phase)           ║"
echo "╠══════════════════════════════════════════════════════════╣"
echo "║  Samples : test/samples_test.tsv                        ║"
echo "║  Phases  : 2 → 5  (Phase 1 skipped — using precomputed) ║"
echo "║  Started : $(date)                    ║"
echo "╚══════════════════════════════════════════════════════════╝"
echo ""

# Check precomputed Phase 1 output exists
CLEANED="results/phase1/stage1_cleaned_sorfs.csv"
SORF_BED="results/phase1/sorfs_genomic.bed"
CDS_BED="results/phase1/cds.bed"

if [ ! -f "${CLEANED}" ]; then
    echo "ERROR: ${CLEANED} not found."
    echo ""
    echo "Download stage1_cleaned_sorfs.csv from the GitHub Releases page"
    echo "and place it at: results/phase1/stage1_cleaned_sorfs.csv"
    echo ""
    echo "Then re-run this script."
    exit 1
fi

# Generate BED files from cleaned sORFs if not already present
mkdir -p results/phase1
if [ ! -f "${SORF_BED}" ]; then
    echo "Generating sORF BED from cleaned table..."
    python3 scripts/make_bed.py \
        --input  "${CLEANED}" \
        --output "${SORF_BED}"
fi

if [ ! -f "${CDS_BED}" ]; then
    echo "Generating CDS BED from GTF..."
    GTF=$(python3 -c "
import yaml
with open('config/config.yaml') as f:
    c = yaml.safe_load(f)
print(c['reference']['gtf'])
")
    python3 scripts/make_cds_bed.py \
        --gtf    "${GTF}" \
        --output "${CDS_BED}"
fi

echo "✓ Phase 1 outputs ready"
echo ""

# Run phases 2–5 with the test sample manifest
bash run_all.sh \
    --samples    test/samples_test.tsv \
    --start_phase 2 \
    --end_phase   5

echo ""
echo "╔══════════════════════════════════════════════════════════╗"
echo "║  TEST COMPLETE                                           ║"
echo "║  Check results/phase5/ for final translated ORF tables  ║"
echo "╚══════════════════════════════════════════════════════════╝"
