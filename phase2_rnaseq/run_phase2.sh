#!/usr/bin/env bash
# =============================================================================
#  run_phase2.sh — RNA-seq alignment (single-end and paired-end)
#
#  Iterates over all RNA-seq samples in samples.tsv and runs:
#    1. SRA download (if SRR provided)
#    2. FastQC (raw)
#    3. Cutadapt trimming
#    4. FastQC (trimmed)
#    5. STAR index build (once)
#    6. STAR alignment → sorted BAM
#
#  Usage:
#    bash phase2_rnaseq/run_phase2.sh [--config config/config.yaml] [--samples config/samples.tsv]
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

GENOME_FA=$(read_cfg reference.genome_fasta)
GTF=$(read_cfg reference.gtf)
STAR_INDEX=$(read_cfg reference.star_index)
THREADS=$(read_cfg threads)
ADAPTER_SE=$(read_cfg rnaseq.adapter_se)
ADAPTER_R1=$(read_cfg rnaseq.adapter_r1)
ADAPTER_R2=$(read_cfg rnaseq.adapter_r2)
QUAL_CUTOFF=$(read_cfg rnaseq.quality_cutoff)
MIN_LEN=$(read_cfg rnaseq.min_length)
OVERHANG_SE=$(read_cfg rnaseq.star_overhang_se)
OVERHANG_PE=$(read_cfg rnaseq.star_overhang_pe)
OUTDIR=$(read_cfg output_dir)

RNA_DATA="data/rna_seq"
mkdir -p "${RNA_DATA}"

PHASE2="${OUTDIR}/phase2"
mkdir -p "${PHASE2}/01_qc/raw" \
         "${PHASE2}/01_qc/trimmed" \
         "${PHASE2}/02_trimmed" \
         "${PHASE2}/03_aligned"

echo "============================================================"
echo "  PHASE 2 — RNA-seq Pipeline"
echo "  Config  : ${CONFIG}"
echo "  Samples : ${SAMPLES}"
echo "  $(date)"
echo "============================================================"

# ── Build STAR index once ─────────────────────────────────────────────────────
if [ ! -d "${STAR_INDEX}" ] || [ -z "$(ls -A ${STAR_INDEX} 2>/dev/null)" ]; then
    echo "── Building STAR index ──"
    mkdir -p "${STAR_INDEX}"
    STAR \
        --runMode genomeGenerate \
        --runThreadN "${THREADS}" \
        --genomeDir "${STAR_INDEX}" \
        --genomeFastaFiles "${GENOME_FA}" \
        --sjdbGTFfile "${GTF}" \
        --sjdbOverhang "${OVERHANG_PE}" \
        --genomeSAindexNbases 14
    echo "✓ STAR index built"
else
    echo "✓ STAR index exists — skipping"
fi

# ── Process each RNA-seq sample ──────────────────────────────────────────────
python3 - <<'PYEOF'
import csv, sys
samples = []
with open("${SAMPLES}") as f:
    for row in csv.DictReader(f, delimiter='\t'):
        if not row['sample_name'].startswith('#') and row['data_type'] == 'rna':
            samples.append(row)
for s in samples:
    print(s['sample_name'])
PYEOF

# Re-read in bash
while IFS=$'\t' read -r sample_name cell_type data_type layout replicate srr fastq_r1 fastq_r2; do
    # Skip header and comments
    [[ "${sample_name}" =~ ^#.*$ || "${sample_name}" == "sample_name" ]] && continue
    [[ "${data_type}" != "rna" ]] && continue

    echo ""
    echo "────────────────────────────────────────────────────────"
    echo "  Sample: ${sample_name}  (${layout})"
    echo "────────────────────────────────────────────────────────"

    # ── Download if SRR provided ─────────────────────────────────────────────
    if [ -n "${srr}" ]; then
        if [ "${layout}" == "single" ]; then
            FQ="${RNA_DATA}/${srr}.fastq.gz"
            if [ ! -f "${FQ}" ]; then
                echo "  Downloading ${srr}..."
                prefetch "${srr}" --output-directory "${RNA_DATA}"
                fasterq-dump "${RNA_DATA}/${srr}" \
                    --outdir "${RNA_DATA}" \
                    --threads "${THREADS}"
                gzip "${RNA_DATA}/${srr}.fastq"
            fi
            fastq_r1="${FQ}"
        else
            FQ_R1="${RNA_DATA}/${srr}_1.fastq.gz"
            FQ_R2="${RNA_DATA}/${srr}_2.fastq.gz"
            if [ ! -f "${FQ_R1}" ]; then
                echo "  Downloading ${srr} (paired-end)..."
                prefetch "${srr}" --output-directory "${RNA_DATA}"
                fasterq-dump "${RNA_DATA}/${srr}" \
                    --outdir "${RNA_DATA}" \
                    --split-files \
                    --threads "${THREADS}"
                gzip "${RNA_DATA}/${srr}_1.fastq"
                gzip "${RNA_DATA}/${srr}_2.fastq"
            fi
            fastq_r1="${FQ_R1}"
            fastq_r2="${FQ_R2}"
        fi
    fi

    SORTED_BAM="${PHASE2}/03_aligned/${sample_name}_sorted.bam"

    # ── FastQC raw ───────────────────────────────────────────────────────────
    if [ ! -f "${PHASE2}/01_qc/raw/${sample_name}_fastqc.zip" ]; then
        echo "  FastQC (raw)..."
        fastqc --outdir "${PHASE2}/01_qc/raw" \
               --threads "${THREADS}" \
               "${fastq_r1}" ${fastq_r2:+"${fastq_r2}"}
    fi

    # ── Cutadapt ─────────────────────────────────────────────────────────────
    if [ "${layout}" == "single" ]; then
        TRIMMED="${PHASE2}/02_trimmed/${sample_name}_trimmed.fastq.gz"
        if [ ! -f "${TRIMMED}" ]; then
            echo "  Trimming (single-end)..."
            cutadapt \
                -a "${ADAPTER_SE}" \
                --quality-cutoff "${QUAL_CUTOFF}" \
                --minimum-length "${MIN_LEN}" \
                --cores "${THREADS}" \
                -o "${TRIMMED}" \
                "${fastq_r1}" \
                2>&1 | tee "${PHASE2}/02_trimmed/${sample_name}_cutadapt.log"
        fi
        TRIM_INPUT="${TRIMMED}"
    else
        TRIMMED_R1="${PHASE2}/02_trimmed/${sample_name}_R1_trimmed.fastq.gz"
        TRIMMED_R2="${PHASE2}/02_trimmed/${sample_name}_R2_trimmed.fastq.gz"
        if [ ! -f "${TRIMMED_R1}" ]; then
            echo "  Trimming (paired-end)..."
            cutadapt \
                -a "${ADAPTER_R1}" -A "${ADAPTER_R2}" \
                --quality-cutoff "${QUAL_CUTOFF}" \
                --minimum-length "${MIN_LEN}" \
                --cores "${THREADS}" \
                -o "${TRIMMED_R1}" -p "${TRIMMED_R2}" \
                "${fastq_r1}" "${fastq_r2}" \
                2>&1 | tee "${PHASE2}/02_trimmed/${sample_name}_cutadapt.log"
        fi
        TRIM_INPUT="${TRIMMED_R1} ${TRIMMED_R2}"
    fi

    # ── FastQC trimmed ───────────────────────────────────────────────────────
    if [ ! -f "${PHASE2}/01_qc/trimmed/${sample_name}_trimmed_fastqc.zip" ]; then
        echo "  FastQC (trimmed)..."
        fastqc --outdir "${PHASE2}/01_qc/trimmed" \
               --threads "${THREADS}" ${TRIM_INPUT}
    fi

    # ── STAR alignment ───────────────────────────────────────────────────────
    if [ ! -f "${SORTED_BAM}" ]; then
        echo "  STAR alignment..."
        STAR \
            --runThreadN "${THREADS}" \
            --genomeDir "${STAR_INDEX}" \
            --readFilesIn ${TRIM_INPUT} \
            --readFilesCommand zcat \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes NH HI AS NM MD \
            --outFileNamePrefix "${PHASE2}/03_aligned/${sample_name}_" \
            --outFilterType BySJout \
            --outFilterMultimapNmax 20 \
            --alignSJoverhangMin 8 \
            --alignSJDBoverhangMin 1 \
            --outFilterMismatchNmax 999 \
            --outFilterMismatchNoverReadLmax 0.04 \
            --alignIntronMin 20 \
            --alignIntronMax 1000000 \
            --alignMatesGapMax 1000000

        mv "${PHASE2}/03_aligned/${sample_name}_Aligned.sortedByCoord.out.bam" \
           "${SORTED_BAM}"
        samtools index "${SORTED_BAM}"
        echo "  ✓ BAM: ${SORTED_BAM}"
    else
        echo "  ✓ BAM exists — skipping alignment"
    fi

done < "${SAMPLES}"

echo ""
echo "============================================================"
echo "  PHASE 2 COMPLETE — $(date)"
echo "============================================================"
