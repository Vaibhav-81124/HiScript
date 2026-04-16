#!/usr/bin/env bash
# =============================================================================
#  run_all.sh — sORF Tool Master Runner
#
#  Runs all phases in order. Edit config/config.yaml and config/samples.tsv
#  before running.
#
#  Usage:
#    bash run_all.sh [OPTIONS]
#
#  Options:
#    --config PATH           Config file (default: config/config.yaml)
#    --samples PATH          Sample manifest (default: config/samples.tsv)
#    --skip_discovery        Use precomputed stage1_novel_sorfs.csv
#    --start_phase N         Start from phase N (1–5). Default: 1
#    --end_phase N           Stop after phase N (1–6). Default: 5 (6 = RF scoring)
# =============================================================================

set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${ROOT}"

CONFIG="config/config.yaml"
SAMPLES="config/samples.tsv"
SKIP_DISC=""
START_PHASE=1
END_PHASE=6

while [[ $# -gt 0 ]]; do
    case $1 in
        --config)         CONFIG="$2";   shift 2 ;;
        --samples)        SAMPLES="$2";  shift 2 ;;
        --skip_discovery) SKIP_DISC="--skip_discovery"; shift ;;
        --start_phase)    START_PHASE="$2"; shift 2 ;;
        --end_phase)      END_PHASE="$2";   shift 2 ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

echo "╔══════════════════════════════════════════════════════════╗"
echo "║              sORF Tool — Full Pipeline                  ║"
echo "╠══════════════════════════════════════════════════════════╣"
echo "║  Config  : ${CONFIG}"
echo "║  Samples : ${SAMPLES}"
echo "║  Phases  : ${START_PHASE} → ${END_PHASE}"
echo "║  Started : $(date)"
echo "╚══════════════════════════════════════════════════════════╝"

run_phase() {
    local PHASE=$1
    local SCRIPT=$2
    shift 2
    if (( PHASE >= START_PHASE && PHASE <= END_PHASE )); then
        echo ""
        echo "▶ Starting Phase ${PHASE}..."
        bash "${SCRIPT}" --config "${CONFIG}" --samples "${SAMPLES}" "$@"
    fi
}

run_phase 1 phase1_discovery/run_phase1.sh  ${SKIP_DISC}
run_phase 2 phase2_rnaseq/run_phase2.sh
run_phase 3 phase3_riboseq/run_phase3.sh
run_phase 4 phase4_translation/run_phase4.sh
run_phase 5 phase5_TE/run_phase5.sh
run_phase 6 phase6_scoring/run_phase6.sh

echo ""
echo "╔══════════════════════════════════════════════════════════╗"
echo "║  ALL PHASES COMPLETE — $(date)"
echo "╚══════════════════════════════════════════════════════════╝"
