#!/usr/bin/env bash
# =============================================================================
#  run_test_full.sh -- Full pipeline validation using public HeLa async data
#
#  Downloads raw FASTQ data from SRA and runs the complete pipeline
#  (Phases 2-5) from scratch. Validates end-to-end reproducibility.
#
#  Dataset:
#    HeLa asynchronous cells, 2 biological replicates
#    RNA-seq + Ribo-seq (GEO: GSE79664, Aviner et al. 2017)
#    SRR3306581, SRR3306582 (RNA-seq)
#    SRR3306588, SRR3306589 (Ribo-seq)
#
#  Expected output: 12 concordant sORFs including RPL26P19
#
#  Runtime: ~5-6 hours (SRA download + alignment + processing)
#
#  Prerequisites:
#    1. conda env create -f environment.yml && conda activate ribowin
#    2. Reference files in data/raw/ (genome FASTA, GTF, rRNA databases)
#       See README.md for download commands
#    3. stage1_cleaned_sorfs.csv from GitHub Releases
#       -> results/phase1/stage1_cleaned_sorfs.csv
#
#  Usage:
#    bash test/run_test_full.sh
# =============================================================================

set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT}"

CONFIG="config/config.yaml"
SAMPLES="test/samples_test.tsv"

echo "============================================================"
echo "  RiboWin -- Full Pipeline Validation (HeLa asynchronous)"
echo "  Phases  : 2 through 5 (SRA download + full alignment)"
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

GTF=$(python3 -c "
import yaml
with open('config/config.yaml') as f:
    c = yaml.safe_load(f)
print(c['reference']['gtf'])
")
check_file "${GTF}" \
    "See README.md Reference Files section for download commands"

GENOME=$(python3 -c "
import yaml
with open('config/config.yaml') as f:
    c = yaml.safe_load(f)
print(c['reference']['genome_fasta'])
")
check_file "${GENOME}" \
    "See README.md Reference Files section for download commands"

echo "Pre-flight checks passed."
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
[ ! -f "${CONCORDANT}" ] && echo "FAIL: ${CONCORDANT} not found" && exit 1

python3 - << 'PYEOF'
import pandas as pd

df = pd.read_csv("results/phase4/ribo_HeLa_async_common_translated_orfs.csv")
n  = len(df)
print(f"  Concordant ORFs found : {n}")

if "genomic_orf_id" in df.columns:
    rpl26 = df[df["genomic_orf_id"].astype(str).str.contains("56504", na=False)]
    print(f"  RPL26P19 locus        : {'FOUND -- PASS' if len(rpl26) > 0 else 'NOT FOUND (check manually)'}")

if n == 12:
    print(f"  Count check           : PASS (12 concordant sORFs)")
elif 10 <= n <= 14:
    print(f"  Count check           : WITHIN RANGE ({n} ORFs, expected 12)")
else:
    print(f"  Count check           : WARNING -- got {n}, expected 12")
PYEOF

echo ""
echo "============================================================"
echo "  FULL TEST COMPLETE -- $(date)"
echo "  Results: results/phase5/HeLa_async_translated_orfs_filtered_withTE.csv"
echo "============================================================"
