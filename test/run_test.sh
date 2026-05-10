#!/usr/bin/env bash
# =============================================================================
#  run_test.sh -- Validation run using precomputed HeLa asynchronous data
#
#  Downloads precomputed files from Zenodo and runs Phases 4+5 only.
#  Runtime: ~10 minutes.
#
#  Prerequisites:
#    1. conda env create -f environment.yml && conda activate ribowin
#    2. GTF in data/raw/ (see README)
#
#  Usage:
#    bash test/run_test.sh
# =============================================================================

set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT}"

CONFIG="config/config.yaml"
SAMPLES="test/samples_test.tsv"
ZENODO_URL="https://zenodo.org/record/20084600/files/RiboWin_test_data.zip"

echo "============================================================"
echo "  RiboWin -- Validation Test (HeLa asynchronous)"
echo "  Started : $(date)"
echo "============================================================"
echo ""

# -- Auto-download Zenodo data if missing -------------------------------------
if [ ! -f "results/phase3/04_psites/HeLa_async_RIBO_rep1_psites.bed" ]; then
    echo "Downloading precomputed test data from Zenodo..."
    wget -q --show-progress "${ZENODO_URL}" -O RiboWin_test_data.zip
    unzip -q RiboWin_test_data.zip
    rm RiboWin_test_data.zip
    echo "Download complete."
    echo ""
fi

# -- Pre-flight checks --------------------------------------------------------
check_file() {
    if [ ! -f "$1" ]; then
        echo "ERROR: Missing required file: $1"
        echo "       $2"
        exit 1
    fi
}

check_file "results/phase1/stage1_cleaned_sorfs.csv" \
    "Should have been downloaded from Zenodo. Re-run this script."

check_file "results/phase2/03_aligned/HeLa_async_RNA_rep1_sorted.bam" \
    "Should have been downloaded from Zenodo. Re-run this script."

check_file "results/phase2/03_aligned/HeLa_async_RNA_rep2_sorted.bam" \
    "Should have been downloaded from Zenodo. Re-run this script."

check_file "results/phase3/04_psites/HeLa_async_RIBO_rep1_psites.bed" \
    "Should have been downloaded from Zenodo. Re-run this script."

check_file "results/phase3/04_psites/HeLa_async_RIBO_rep2_psites.bed" \
    "Should have been downloaded from Zenodo. Re-run this script."

GTF=$(python3 -c "
import yaml
with open('config/config.yaml') as f:
    c = yaml.safe_load(f)
print(c['reference']['gtf'])
")

check_file "${GTF}" \
    "Download GTF: wget -P data/raw/ https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz && gunzip data/raw/Homo_sapiens.GRCh38.115.gtf.gz"

echo "All required files present."
echo ""

# -- Generate BED files -------------------------------------------------------
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

# -- Verify chr consistency before running ------------------------------------
echo "Checking chromosome naming consistency..."

BED_CHR=$(head -1 results/phase1/sorfs_genomic.bed | cut -f1)
PSITE_CHR=$(head -1 results/phase3/04_psites/HeLa_async_RIBO_rep1_psites.bed | cut -f1)

echo "  sORF BED chrom:  ${BED_CHR}"
echo "  P-site BED chrom: ${PSITE_CHR}"

if [[ "${BED_CHR}" == chr* ]] && [[ "${PSITE_CHR}" != chr* ]]; then
    echo "  WARNING: chr prefix mismatch detected. Stripping chr from sORF BED..."
    sed -i 's/^chr//' results/phase1/sorfs_genomic.bed
    sed -i 's/^chr//' results/phase1/cds.bed
    echo "  Fixed."
elif [[ "${BED_CHR}" != chr* ]] && [[ "${PSITE_CHR}" == chr* ]]; then
    echo "  WARNING: chr prefix mismatch detected. Stripping chr from P-site BEDs..."
    sed -i 's/^chr//' results/phase3/04_psites/HeLa_async_RIBO_rep1_psites.bed
    sed -i 's/^chr//' results/phase3/04_psites/HeLa_async_RIBO_rep2_psites.bed
    echo "  Fixed."
else
    echo "  OK: chromosome naming is consistent."
fi
echo ""

# -- Run phases 4 and 5 -------------------------------------------------------
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
import pandas as pd

df = pd.read_csv("results/phase4/ribo_HeLa_async_common_translated_orfs.csv")
n  = len(df)
print(f"  Concordant ORFs found : {n}")

if "genomic_orf_id" in df.columns:
    rpl26 = df[df["genomic_orf_id"].astype(str).str.contains("56504", na=False)]
    if len(rpl26) > 0:
        print(f"  RPL26P19 locus        : FOUND -- PASS")
    else:
        print(f"  RPL26P19 locus        : not found by coordinate (check manually)")

if n == 12:
    print(f"  Count check           : PASS (12 concordant ORFs)")
elif 8 <= n <= 16:
    print(f"  Count check           : WITHIN RANGE ({n} ORFs, expected ~12)")
else:
    print(f"  Count check           : WARNING -- got {n}, expected ~12")
PYEOF

echo ""
echo "============================================================"
echo "  TEST COMPLETE -- $(date)"
echo "  Results: results/phase5/HeLa_async_translated_orfs_filtered_withTE.csv"
echo "============================================================"
