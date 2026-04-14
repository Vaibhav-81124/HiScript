#!/usr/bin/env bash
# =============================================================================
#  run_phase5.sh — Translation Efficiency (per cell type)
#
#  For each cell type:
#    1. TE_compute.py  — merge ribo + RNA counts, compute TE with pseudocount
#    2. TE_filter.py   — restrict to translated ORFs, recalculate TE cleanly
#
#  Usage:
#    bash phase5_TE/run_phase5.sh [--config config/config.yaml] [--samples config/samples.tsv]
# =============================================================================

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${ROOT}"

CONFIG="config/config.yaml"
SAMPLES="config/samples.tsv"

while [[ $# -gt 0 ]]; do
    case $1 in
        --config)  CONFIG="$2";  shift 2 ;;
        --samples) SAMPLES="$2"; shift 2 ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

read_cfg() {
    python3 -c "
import yaml
with open('${CONFIG}') as f:
    c = yaml.safe_load(f)
keys = '$1'.split('.')
v = c
for k in keys:
    v = v[k]
print(v)
"
}

PSEUDO=$(read_cfg te.pseudocount)
OUTDIR=$(read_cfg output_dir)

PHASE4="${OUTDIR}/phase4"
PHASE5="${OUTDIR}/phase5"
mkdir -p "${PHASE5}"

echo "============================================================"
echo "  PHASE 5 — Translation Efficiency"
echo "  Config  : ${CONFIG}"
echo "  Samples : ${SAMPLES}"
echo "  $(date)"
echo "============================================================"

mapfile -t CELL_TYPES < <(
    awk -F'\t' 'NR>1 && $1!~/^#/ {print $2}' "${SAMPLES}" | sort -u
)

for CELL in "${CELL_TYPES[@]}"; do

    echo ""
    echo "────────────────────────────────────────────────────────"
    echo "  Cell type: ${CELL}"
    echo "────────────────────────────────────────────────────────"

    # Resolve per-rep count files
    RIBO1="${PHASE4}/${CELL}_rep1_ribo_counts.txt"
    RIBO2="${PHASE4}/${CELL}_rep2_ribo_counts.txt"
    RNA1="${PHASE4}/${CELL}_rep1_rna_counts_clean.txt"
    RNA2="${PHASE4}/${CELL}_rep2_rna_counts_clean.txt"
    MERGED="${PHASE4}/ribo_${CELL}_common_translated_orfs.csv"

    for F in "${RIBO1}" "${RIBO2}" "${RNA1}" "${RNA2}" "${MERGED}"; do
        if [ ! -f "${F}" ]; then
            echo "  ERROR: Missing required file: ${F}"
            echo "  Run phase 4 first."
            continue 2
        fi
    done

    # ── TE_compute ────────────────────────────────────────────────────────────
    TE_TABLE="${PHASE5}/${CELL}_translation_efficiency.csv"
    if [ ! -f "${TE_TABLE}" ]; then
        echo "  Computing TE..."
        python3 scripts/TE_compute.py \
            --ribo1   "${RIBO1}" \
            --ribo2   "${RIBO2}" \
            --rna1    "${RNA1}" \
            --rna2    "${RNA2}" \
            --cell    "${CELL}" \
            --pseudo  "${PSEUDO}" \
            --outdir  "${PHASE5}"
        echo "  ✓ TE table: ${TE_TABLE}"
    fi

    # ── TE_filter ─────────────────────────────────────────────────────────────
    FILTERED="${PHASE5}/${CELL}_translated_orfs_filtered_withTE.csv"
    if [ ! -f "${FILTERED}" ]; then
        echo "  Filtering to translated ORFs..."
        python3 scripts/TE_filter.py \
            --te_table   "${TE_TABLE}" \
            --translated "${MERGED}" \
            --cell       "${CELL}" \
            --outdir     "${PHASE5}"
        echo "  ✓ Filtered: ${FILTERED}"
    fi

done

echo ""
echo "============================================================"
echo "  PHASE 5 COMPLETE — $(date)"
echo "  Final translated ORFs per cell type in: ${PHASE5}/"
echo "============================================================"
