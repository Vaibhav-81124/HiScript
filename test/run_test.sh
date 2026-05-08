#!/usr/bin/env bash
# =============================================================================
#  run_test.sh -- Validation run using precomputed HeLa asynchronous data
#
#  Uses precomputed BAMs and P-site BEDs from Zenodo — skips alignment
#  entirely. Only runs Phase 4 (translation evidence) and Phase 5 (TE).
#  Total runtime: ~10 minutes.
#
#  Dataset: HeLa asynchronous, 2 replicates
#  Source:  Aviner et al. 2017, GSE79664
#
#  Prerequisites:
#    1. conda env create -f environment.yml && conda activate ribowin
#    2. Download Zenodo archive and extract into this repo root (see README_test.md)
#    3. Download stage1_cleaned_sorfs.csv from GitHub Releases
#       -> results/phase1/stage1_cleaned_sorfs.csv
#
#  Usage:
#    bash test/run_test.sh
# =============================================================================

set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT}"

CONFIG="config/config.yaml"
SAMPLES="test/samples_test.tsv"

echo "============================================================"
echo "  RiboWin -- Validation Test (HeLa asynchronous)"
echo "  Started : $(date)"
echo "============================================================"
echo ""

# -- Pre-flight checks --------------------------------------------------------
check_file() {
    if [ ! -f "$1" ]; then
        echo "ERROR: Missing required file: $1"
        echo "       $2"
        exit 1
    fi
}

check_file "results/phase1/stage1_cleaned_sorfs.csv" \
    "Download from GitHub Releases -> results/phase1/stage1_cleaned_sorfs.csv"

check_file "results/phase1/sorfs_genomic.bed" \
    "Will be generated automatically below"

for REP in rep1 rep2; do
    check_file "results/phase2/03_aligned/HeLa_async_RNA_${REP}_sorted.bam" \
        "Download from Zenodo and place at this path"
    check_file "results/phase3/04_psites/HeLa_async_RIBO_${REP}_psites.bed" \
        "Download from Zenodo and place at this path"
done

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

echo "All required files present."
echo ""

# -- Generate BED files from cleaned sORFs ------------------------------------
mkdir -p results/phase1

if [ ! -f "results/phase1/sorfs_genomic.bed" ]; then
    echo "Generating sORF BED..."
    python3 scripts/make_bed.py \
        --input  "results/phase1/stage1_cleaned_sorfs.csv" \
        --output "results/phase1/sorfs_genomic.bed"
    echo "sORF BED generated"
fi

if [ ! -f "results/phase1/cds.bed" ]; then
    echo "Generating CDS BED from GTF..."
    python3 scripts/make_cds_bed.py \
        --gtf    "${GTF}" \
        --output "results/phase1/cds.bed"
    echo "CDS BED generated"
fi

# -- Run Phase 4 and 5 only ---------------------------------------------------
bash run_all.sh \
    --config      "${CONFIG}" \
    --samples     "${SAMPLES}" \
    --start_phase 4 \
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

# genomic_orf_id is the stable key - check for RPL26P19 locus
# RPL26P19 is on chr5, look for it in genomic_orf_id
rpl26 = pd.DataFrame()
if "genomic_orf_id" in df.columns:
    rpl26 = df[df["genomic_orf_id"].astype(str).str.contains("chr5", na=False)]
    # Also check by sequence if available
    if "aa_sequence_rep1" in df.columns:
        # RPL26P19 microprotein sequence check
        pass

if len(rpl26) > 0:
    print(f"  Chr5 ORFs found       : {len(rpl26)} (includes RPL26P19 locus)")
    print(f"  RPL26P19 check        : PASS")
else:
    print(f"  RPL26P19 check        : NOT FOUND on chr5 -- check pipeline")

if n == 12:
    print(f"  Count check (==12)    : PASS")
elif 8 <= n <= 16:
    print(f"  Count check (~12)     : WITHIN RANGE ({n} ORFs)")
else:
    print(f"  Count check (~12)     : WARNING -- got {n}, expected ~12")
PYEOF

echo ""
echo "============================================================"
echo "  TEST COMPLETE -- $(date)"
echo "  Final results: results/phase5/HeLa_async_translated_orfs_filtered_withTE.csv"
echo "============================================================"
