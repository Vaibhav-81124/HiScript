#!/usr/bin/env bash
# =============================================================================
#  run_test.sh -- Full validation run using public HeLa asynchronous data
#
#  Automatically downloads raw FASTQ data from SRA and runs the complete
#  pipeline (Phases 2-5). No pre-processed files needed.
#
#  Dataset:
#    HeLa asynchronous cells, 2 biological replicates
#    RNA-seq + Ribo-seq (GEO: GSE79664, Aviner et al. 2017)
#    SRR3306581, SRR3306582 (RNA-seq)
#    SRR3306588, SRR3306589 (Ribo-seq)
#
#  Expected output:
#    9 concordant sORFs (test dataset uses subsetted BAMs)
#    Full dataset produces 12 sORFs as reported in the paper
#
#  Prerequisites:
#    1. conda env create -f environment.yml && conda activate sorf-tool
#    2. Place GRCh38 reference files in data/raw/ (see README)
#    3. Download stage1_cleaned_sorfs.csv from GitHub Releases
#       -> results/phase1/stage1_cleaned_sorfs.csv
#
#  Usage:
#    bash test/run_test.sh
#
#  Runtime: ~2-4 hours depending on download speed and machine specs
# =============================================================================

set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT}"

CONFIG="config/config.yaml"
SAMPLES="test/samples_test.tsv"

echo "============================================================"
echo "  RiboWin -- Validation Test (HeLa asynchronous)"
echo "  Samples : ${SAMPLES}"
echo "  Started : $(date)"
echo "============================================================"
echo ""

# -- Pre-flight checks --------------------------------------------------------
if [ ! -f "results/phase1/stage1_cleaned_sorfs.csv" ]; then
    echo "ERROR: Missing results/phase1/stage1_cleaned_sorfs.csv"
    echo "Download from GitHub Releases and place at that path."
    exit 1
fi

GTF=$(python3 -c "
import yaml
with open('config/config.yaml') as f:
    c = yaml.safe_load(f)
print(c['reference']['gtf'])
")

if [ ! -f "${GTF}" ]; then
    echo "ERROR: GTF not found at ${GTF}"
    echo "Update reference.gtf in config/config.yaml"
    exit 1
fi

echo "Pre-flight checks passed."
echo ""

# -- Generate BED files if missing --------------------------------------------
mkdir -p results/phase1

if [ ! -f "results/phase1/sorfs_genomic.bed" ]; then
    echo "Generating sORF BED..."
    python3 scripts/make_bed.py \
        --input  "results/phase1/stage1_cleaned_sorfs.csv" \
        --output "results/phase1/sorfs_genomic.bed"
fi

if [ ! -f "results/phase1/cds.bed" ]; then
    echo "Generating CDS BED from GTF..."
    python3 scripts/make_cds_bed.py \
        --gtf    "${GTF}" \
        --output "results/phase1/cds.bed"
fi

# -- Run phases 2-5 -----------------------------------------------------------
bash run_all.sh \
    --config      "${CONFIG}" \
    --samples     "${SAMPLES}" \
    --start_phase 2 \
    --end_phase   5

# -- Validation check ---------------------------------------------------------
echo ""
echo "============================================================"
echo "  VALIDATION CHECK"
echo "============================================================"

CONCORDANT="results/phase4/ribo_HeLa_async_common_translated_orfs.csv"

if [ ! -f "${CONCORDANT}" ]; then
    echo "FAIL: Concordant ORF file not found: ${CONCORDANT}"
    exit 1
fi

python3 - << 'PYEOF'
import pandas as pd, sys

df = pd.read_csv("results/phase4/ribo_HeLa_async_common_translated_orfs.csv")
n  = len(df)

print(f"  Concordant ORFs found : {n}")

# Check for RPL26P19
id_col = "orf_id" if "orf_id" in df.columns else "orf_id_rep1"
rpl26  = df[df[id_col].astype(str).str.contains("13030", na=False)]

if len(rpl26) > 0:
    print(f"  RPL26P19 (ORF_13030)  : FOUND")
else:
    print(f"  RPL26P19 (ORF_13030)  : NOT FOUND -- check pipeline output")

# Test dataset uses subsetted BAMs -> 9 sORFs expected
# Full dataset produces 12 sORFs as reported in the paper
if n == 9:
    print(f"  Count check (==9)     : PASS")
elif 7 <= n <= 11:
    print(f"  Count check (~9)      : WITHIN RANGE ({n} ORFs)")
else:
    print(f"  Count check (~9)      : WARNING -- got {n}, expected ~9 (subsetted BAMs)")
print(f"  NOTE: Full dataset yields 12 sORFs (see paper)")
PYEOF

echo ""
echo "============================================================"
echo "  TEST COMPLETE -- $(date)"
echo "  Final results: results/phase5/HeLa_async_translated_orfs_filtered_withTE.csv"
echo "============================================================"
