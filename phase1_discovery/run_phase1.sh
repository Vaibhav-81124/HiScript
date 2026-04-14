#!/usr/bin/env bash
# =============================================================================
#  run_phase1.sh — sORF Discovery
#
#  Runs:
#    1a. sorf_discovery.py   → raw sORF table
#    1b. sorf_cleaning.py    → deduplicated + filtered sORF table
#    1c. make_bed.py         → genomic BED file
#    1d. make_cds_bed.py     → CDS BED (used in Phase 4)
#
#  Usage:
#    bash phase1_discovery/run_phase1.sh [--config config/config.yaml] [--skip_discovery]
#
#  Flags:
#    --config PATH           Path to config.yaml (default: config/config.yaml)
#    --skip_discovery        Skip step 1a and use precomputed stage1_novel_sorfs.csv
# =============================================================================

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${ROOT}"

# ── Defaults ────────────────────────────────────────────────────────────────
CONFIG="config/config.yaml"
SKIP_DISCOVERY=false

# ── Argument parsing ─────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case $1 in
        --config)           CONFIG="$2";        shift 2 ;;
        --skip_discovery)   SKIP_DISCOVERY=true; shift ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

# ── Read config via Python (avoids yq dependency) ────────────────────────────
read_cfg() {
    python3 -c "
import yaml, sys
with open('${CONFIG}') as f:
    c = yaml.safe_load(f)
keys = '$1'.split('.')
v = c
for k in keys:
    v = v[k]
print(v)
"
}

CDNA_FASTA=$(read_cfg reference.cdna_fasta)
GTF=$(read_cfg reference.gtf)
MIN_AA=$(read_cfg discovery.min_aa)
MAX_AA=$(read_cfg discovery.max_aa)
MIN_AA_CLEAN=$(read_cfg discovery.min_aa_after_cleaning)
OUTDIR=$(read_cfg output_dir)

PHASE1_OUT="${OUTDIR}/phase1"
mkdir -p "${PHASE1_OUT}"

echo "============================================================"
echo "  PHASE 1 — sORF Discovery"
echo "  Config : ${CONFIG}"
echo "  Output : ${PHASE1_OUT}"
echo "  $(date)"
echo "============================================================"

# ── Step 1a: Discovery ───────────────────────────────────────────────────────
RAW_CSV="${PHASE1_OUT}/stage1_novel_sorfs.csv"

if [ "${SKIP_DISCOVERY}" = true ]; then
    echo "⏭  Skipping discovery (--skip_discovery set)"
    if [ ! -f "${RAW_CSV}" ]; then
        echo "ERROR: --skip_discovery requires ${RAW_CSV} to exist"
        exit 1
    fi
else
    echo "── Step 1a: sORF scanning ──"
    python3 scripts/sorf_discovery.py \
        --fasta   "${CDNA_FASTA}" \
        --gtf     "${GTF}" \
        --min_aa  "${MIN_AA}" \
        --max_aa  "${MAX_AA}" \
        --output  "${RAW_CSV}"
    echo "✓ Discovery complete"
fi

# ── Step 1b: Cleaning ────────────────────────────────────────────────────────
CLEAN_CSV="${PHASE1_OUT}/stage1_cleaned_sorfs.csv"

echo "── Step 1b: Cleaning and deduplication ──"
python3 scripts/sorf_cleaning.py \
    --input   "${RAW_CSV}" \
    --min_aa  "${MIN_AA_CLEAN}" \
    --output  "${CLEAN_CSV}"
echo "✓ Cleaning complete"

# ── Step 1c: BED conversion ──────────────────────────────────────────────────
SORF_BED="${PHASE1_OUT}/sorfs_genomic.bed"

echo "── Step 1c: BED file generation ──"
python3 scripts/make_bed.py \
    --input   "${CLEAN_CSV}" \
    --output  "${SORF_BED}"
echo "✓ BED file written"

# ── Step 1d: CDS BED ─────────────────────────────────────────────────────────
CDS_BED="${PHASE1_OUT}/cds.bed"

echo "── Step 1d: CDS BED from GTF ──"
python3 scripts/make_cds_bed.py \
    --gtf     "${GTF}" \
    --output  "${CDS_BED}"
echo "✓ CDS BED written"

echo "============================================================"
echo "  PHASE 1 COMPLETE"
echo "  Cleaned sORFs : ${CLEAN_CSV}"
echo "  sORF BED      : ${SORF_BED}"
echo "  CDS BED       : ${CDS_BED}"
echo "============================================================"
