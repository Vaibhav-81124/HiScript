#!/usr/bin/env bash
# =============================================================================
#  run_phase6.sh -- RF scoring of translated sORFs
#
#  For each cell type, takes the Phase 5 TE-filtered table and runs:
#    orf_scorer.py  -> rf_score (0-1) + rf_label per ORF
#
#  First cell type trains the model (self-supervised).
#  Subsequent cell types can reuse it with --reuse_model, or train their own.
#
#  Usage:
#    bash phase6_scoring/run_phase6.sh [--config config/config.yaml]
#                                      [--samples config/samples.tsv]
#                                      [--reuse_model]
# =============================================================================

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${ROOT}"

CONFIG="config/config.yaml"
SAMPLES="config/samples.tsv"
REUSE_MODEL=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --config)       CONFIG="$2";       shift 2 ;;
        --samples)      SAMPLES="$2";      shift 2 ;;
        --reuse_model)  REUSE_MODEL=true;  shift ;;
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

OUTDIR=$(read_cfg output_dir)
PHASE5="${OUTDIR}/phase5"
PHASE6="${OUTDIR}/phase6_scoring"
mkdir -p "${PHASE6}"

echo "============================================================"
echo "  PHASE 6 -- RF ORF Scoring"
echo "  Config  : ${CONFIG}"
echo "  Samples : ${SAMPLES}"
echo "  $(date)"
echo "============================================================"

mapfile -t CELL_TYPES < <(
    awk -F'\t' 'NR>1 && $1!~/^#/ {print $2}' "${SAMPLES}" | sort -u
)

FIRST_MODEL=""

for CELL in "${CELL_TYPES[@]}"; do

    echo ""
    echo "------------------------------------------------------------"
    echo "  Cell type: ${CELL}"
    echo "------------------------------------------------------------"

    INPUT="${PHASE5}/${CELL}_translated_orfs_filtered_withTE.csv"

    if [ ! -f "${INPUT}" ]; then
        echo "  WARNING: No Phase 5 output found for ${CELL} -- skipping"
        echo "  Expected: ${INPUT}"
        continue
    fi

    SCORED="${PHASE6}/${CELL}_scored_orfs.csv"

    if [ ! -f "${SCORED}" ]; then

        MODEL_ARG=""
        if [ "${REUSE_MODEL}" = true ] && [ -n "${FIRST_MODEL}" ]; then
            echo "  Using pre-trained model from: ${FIRST_MODEL}"
            MODEL_ARG="--model ${FIRST_MODEL}"
        fi

        python3 scripts/orf_scorer.py \
            --input   "${INPUT}" \
            --cell    "${CELL}" \
            --outdir  "${PHASE6}" \
            ${MODEL_ARG}

        echo "  Scored: ${SCORED}"
    else
        echo "  Scored file exists -- skipping"
    fi

    # Track first trained model for reuse
    if [ -z "${FIRST_MODEL}" ]; then
        FIRST_MODEL="${PHASE6}/${CELL}_rf_model.pkl"
    fi

done

echo ""
echo "============================================================"
echo "  PHASE 6 COMPLETE -- $(date)"
echo "  Scored tables in: ${PHASE6}/"
echo "  Key output columns added: rf_score, rf_label"
echo "============================================================"
