#!/usr/bin/env bash
# =============================================================================
#  run_phase4.sh — Translation Evidence (per cell type)
#
#  For each cell type in samples.tsv, runs across both replicates:
#    1.  BED → GTF for featureCounts
#    2.  featureCounts  (RNA counts per sORF)
#    3.  bedtools intersect (Ribo P-sites × sORF BED)
#    4.  awk aggregate  (Ribo counts per ORF)
#    5.  stage4_translation.py  (CDS removal + read filter)
#    6.  stage5_periodicity.py  (triplet periodicity filter)
#    7.  merge_reps.py          (keep reproducible ORFs across reps)
#
#  Usage:
#    bash phase4_translation/run_phase4.sh [--config config/config.yaml] [--samples config/samples.tsv]
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

sed -i 's/\r$//' "${SAMPLES}" 2>/dev/null || true

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

THREADS=$(read_cfg threads)
MIN_READS=$(read_cfg translation.min_ribo_reads)
MICRO_MAX=$(read_cfg translation.micro_max_aa)
MIN_PSITE=$(read_cfg translation.min_psite_reads)
FRAME_THR=$(read_cfg translation.frame0_threshold)
OUTDIR=$(read_cfg output_dir)

PHASE2="${OUTDIR}/phase2"
PHASE3="${OUTDIR}/phase3"
PHASE4="${OUTDIR}/phase4"

SORF_BED="results/phase1/sorfs_genomic.bed"
CDS_BED="results/phase1/cds.bed"
SORFS_CSV="results/phase1/stage1_cleaned_sorfs.csv"

mkdir -p "${PHASE4}"

# ── GTF for featureCounts (derived once from sORF BED) ───────────────────────
SORF_GTF="${PHASE4}/sorfs_annotation.gtf"
if [ ! -f "${SORF_GTF}" ]; then
    echo "── Generating sORF annotation GTF ──"
    awk 'BEGIN{OFS="\t"} {
        if ($2 < $3)
            print $1,"sorf_project","exon",$2+1,$3,".",$6,".",\
                  "gene_id \""$4"\"; transcript_id \""$4"\";"
    }' "${SORF_BED}" > "${SORF_GTF}"
    echo "✓ GTF: ${SORF_GTF}"
fi

echo "============================================================"
echo "  PHASE 4 — Translation Evidence"
echo "  Config  : ${CONFIG}"
echo "  Samples : ${SAMPLES}"
echo "  $(date)"
echo "============================================================"

# ── Collect unique cell types ─────────────────────────────────────────────────
mapfile -t CELL_TYPES < <(
    awk -F'\t' 'NR>1 && $1!~/^#/ {print $2}' "${SAMPLES}" | sort -u
)

for CELL in "${CELL_TYPES[@]}"; do

    echo ""
    echo "────────────────────────────────────────────────────────"
    echo "  Cell type: ${CELL}"
    echo "────────────────────────────────────────────────────────"

    # Collect ribo samples and rna samples for this cell type
    mapfile -t RIBO_SAMPLES < <(
        awk -F'\t' -v cell="${CELL}" 'NR>1 && $1!~/^#/ && $2==cell && $3=="ribo" {print $1"__"$5}' "${SAMPLES}"
    )
    mapfile -t RNA_SAMPLES < <(
        awk -F'\t' -v cell="${CELL}" 'NR>1 && $1!~/^#/ && $2==cell && $3=="rna" {print $1"__"$5}' "${SAMPLES}"
    )

    # Per-replicate processing
    STAGE5_FILES=()

    for RIBO_ENTRY in "${RIBO_SAMPLES[@]}"; do
        RIBO_SAMPLE="${RIBO_ENTRY%%__*}"
        REP="${RIBO_ENTRY##*__}"

        echo "  Processing ribo replicate: ${RIBO_SAMPLE} (${REP})"

        PSITE_BED="${PHASE3}/04_psites/${RIBO_SAMPLE}_psites.bed"
        RIBO_SEG="${PHASE4}/${CELL}_${REP}_segment_counts.bed"
        RIBO_TOT="${PHASE4}/${CELL}_${REP}_ribo_counts.txt"

        # ── Ribo counts via bedtools ──────────────────────────────────────────
        if [ ! -f "${RIBO_TOT}" ]; then
            echo "    bedtools intersect..."
            bedtools intersect \
                -a "${SORF_BED}" \
                -b "${PSITE_BED}" \
                -s -c > "${RIBO_SEG}"

            awk '{counts[$4]+=$7} END {for (orf in counts) print orf"\t"counts[orf]}' \
                "${RIBO_SEG}" > "${RIBO_TOT}"
            echo "    ✓ Ribo counts: ${RIBO_TOT}"
        fi

        # ── RNA counts for matching RNA rep ──────────────────────────────────
        # Find RNA sample for same replicate
        RNA_MATCH=""
        for RNA_ENTRY in "${RNA_SAMPLES[@]}"; do
            RNA_SAMPLE="${RNA_ENTRY%%__*}"
            RNA_REP="${RNA_ENTRY##*__}"
            if [[ "${RNA_REP}" == "${REP}" ]]; then
                RNA_MATCH="${RNA_SAMPLE}"
                break
            fi
        done

        RNA_BAM="${PHASE2}/03_aligned/${RNA_MATCH}_sorted.bam"

        # Detect paired-end from samples TSV
        IS_PAIRED=$(awk -F'\t' -v s="${RNA_MATCH}" '$1==s {print $4}' "${SAMPLES}")

        RNA_RAW="${PHASE4}/${CELL}_${REP}_rna_counts_raw.txt"
        RNA_CLEAN="${PHASE4}/${CELL}_${REP}_rna_counts_clean.txt"

        if [ ! -f "${RNA_CLEAN}" ] && [ -n "${RNA_MATCH}" ]; then
            echo "    featureCounts (RNA)..."
            if [[ "${IS_PAIRED}" == "paired" ]]; then
                featureCounts \
                    -T "${THREADS}" \
                    -a "${SORF_GTF}" \
                    -o "${RNA_RAW}" \
                    -t exon -g gene_id \
                    -O -s 0 -p --countReadPairs \
                    "${RNA_BAM}"
            else
                featureCounts \
                    -T "${THREADS}" \
                    -a "${SORF_GTF}" \
                    -o "${RNA_RAW}" \
                    -t exon -g gene_id \
                    "${RNA_BAM}"
            fi
            awk 'NR>2 {print $1"\t"$7}' "${RNA_RAW}" > "${RNA_CLEAN}"
            echo "    ✓ RNA counts: ${RNA_CLEAN}"
        fi

        # ── Stage 4 translation filter ────────────────────────────────────────
        S4_TAG="${CELL}_${REP}"
        S4_HCO="${PHASE4}/stage4_${S4_TAG}_high_confidence_novel_orfs.csv"

        if [ ! -f "${S4_HCO}" ]; then
            echo "    Stage 4 translation filter..."
            python3 scripts/stage4_translation.py \
                --sorfs     "${SORFS_CSV}" \
                --ribo      "${RIBO_TOT}" \
                --cds_bed   "${CDS_BED}" \
                --sorf_bed  "${SORF_BED}" \
                --sample    "${S4_TAG}" \
                --min_reads "${MIN_READS}" \
                --micro_max "${MICRO_MAX}" \
                --outdir    "${PHASE4}"
            echo "    ✓ Stage 4 complete"
        fi

        # ── Stage 5 periodicity filter ────────────────────────────────────────
        S5_HC="${PHASE4}/stage5_${S4_TAG}_high_confidence_translated_orfs.csv"

        if [ ! -f "${S5_HC}" ]; then
            echo "    Stage 5 periodicity filter..."
            python3 scripts/stage5_periodicity.py \
                --orfs      "${S4_HCO}" \
                --psites    "${PSITE_BED}" \
                --sample    "${S4_TAG}" \
                --min_reads "${MIN_PSITE}" \
                --frame_thr "${FRAME_THR}" \
                --outdir    "${PHASE4}"
            echo "    ✓ Stage 5 complete"
        fi

        STAGE5_FILES+=("${S5_HC}")

    done

    # ── Merge replicates ──────────────────────────────────────────────────────
    MERGED="${PHASE4}/ribo_${CELL}_common_translated_orfs.csv"
    if [ ! -f "${MERGED}" ] && [ "${#STAGE5_FILES[@]}" -ge 2 ]; then
        echo "  Merging replicates for ${CELL}..."
        python3 scripts/merge_reps.py \
            --rep1   "${STAGE5_FILES[0]}" \
            --rep2   "${STAGE5_FILES[1]}" \
            --cell   "${CELL}" \
            --outdir "${PHASE4}"
        echo "  ✓ Merged: ${MERGED}"
    fi

done

echo ""
echo "============================================================"
echo "  PHASE 4 COMPLETE — $(date)"
echo "============================================================"
