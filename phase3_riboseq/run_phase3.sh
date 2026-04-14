#!/usr/bin/env bash
# =============================================================================
#  run_phase3.sh — Ribo-seq processing
#
#  Runs for each ribo sample in samples.tsv:
#    1. SRA download (if SRR provided)
#    2. Cutadapt trim + size-select (RPF_MIN–RPF_MAX)
#    3. rRNA removal (SortMeRNA)
#    4. STAR alignment
#    5. RPF length filter (samtools)
#    6. P-site assignment → BED
#
#  Usage:
#    bash phase3_riboseq/run_phase3.sh [--config config/config.yaml] [--samples config/samples.tsv]
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

GTF=$(read_cfg reference.gtf)
STAR_INDEX=$(read_cfg reference.star_index)
RRNA_DB=$(read_cfg reference.rrna_db)
THREADS=$(read_cfg threads)
ADAPTER=$(read_cfg riboseq.adapter)
RPF_MIN=$(read_cfg riboseq.rpf_min)
RPF_MAX=$(read_cfg riboseq.rpf_max)
QUAL_CUTOFF=$(read_cfg riboseq.quality_cutoff)
OUTDIR=$(read_cfg output_dir)

RIBO_DATA="data/ribo_seq"
mkdir -p "${RIBO_DATA}"

PHASE3="${OUTDIR}/phase3"
mkdir -p "${PHASE3}/01_trimmed" \
         "${PHASE3}/02_norrna" \
         "${PHASE3}/03_aligned" \
         "${PHASE3}/04_psites"

echo "============================================================"
echo "  PHASE 3 — Ribo-seq Pipeline"
echo "  Config  : ${CONFIG}"
echo "  Samples : ${SAMPLES}"
echo "  $(date)"
echo "============================================================"

# ── Validate STAR index ───────────────────────────────────────────────────────
if [ ! -f "${STAR_INDEX}/SA" ]; then
    echo "ERROR: STAR index not found at ${STAR_INDEX}. Run phase 2 first."
    exit 1
fi

while IFS=$'\t' read -r sample_name cell_type data_type layout replicate srr fastq_r1 fastq_r2; do
    [[ "${sample_name}" =~ ^#.*$ || "${sample_name}" == "sample_name" ]] && continue
    [[ "${data_type}" != "ribo" ]] && continue

    echo ""
    echo "────────────────────────────────────────────────────────"
    echo "  Sample: ${sample_name}"
    echo "────────────────────────────────────────────────────────"

    # ── Download ─────────────────────────────────────────────────────────────
    if [ -n "${srr}" ]; then
        FQ="${RIBO_DATA}/${srr}.fastq.gz"
        if [ ! -f "${FQ}" ]; then
            echo "  Downloading ${srr}..."
            prefetch "${srr}" --output-directory "${RIBO_DATA}"
            fasterq-dump "${RIBO_DATA}/${srr}" \
                --outdir "${RIBO_DATA}" \
                --threads "${THREADS}"
            gzip "${RIBO_DATA}/${srr}.fastq"
        fi
        fastq_r1="${FQ}"
    fi

    TRIMMED="${PHASE3}/01_trimmed/${sample_name}_trimmed.fastq.gz"
    NORRNA="${PHASE3}/02_norrna/${sample_name}_norRNA.fq.gz"
    RPF_BAM="${PHASE3}/03_aligned/${sample_name}_RPFs.bam"
    PSITE_BED="${PHASE3}/04_psites/${sample_name}_psites.bed"

    # ── Step 1: Trim + size select ───────────────────────────────────────────
    if [ ! -f "${TRIMMED}" ]; then
        echo "  Trimming and size-selecting (${RPF_MIN}–${RPF_MAX} nt)..."
        cutadapt \
            -a "${ADAPTER}" \
            --quality-cutoff "${QUAL_CUTOFF}" \
            --minimum-length "${RPF_MIN}" \
            --maximum-length "${RPF_MAX}" \
            --cores "${THREADS}" \
            -o "${TRIMMED}" \
            "${fastq_r1}"
        echo "  ✓ Trimmed"
    fi

    # ── Step 2: rRNA removal ─────────────────────────────────────────────────
    if [ ! -f "${NORRNA}" ]; then
        echo "  rRNA removal (SortMeRNA)..."
        WORKDIR="${PHASE3}/02_norrna/sortmerna_${sample_name}_workdir"
        rm -rf "${WORKDIR}" && mkdir -p "${WORKDIR}"

        sortmerna \
            --ref "${RRNA_DB}/silva-euk-28s-id98.fasta" \
            --ref "${RRNA_DB}/silva-euk-18s-id95.fasta" \
            --ref "${RRNA_DB}/rfam-5s-database-id98.fasta" \
            --ref "${RRNA_DB}/rfam-5.8s-database-id98.fasta" \
            --reads "${TRIMMED}" \
            --workdir "${WORKDIR}" \
            --other "${PHASE3}/02_norrna/${sample_name}_norRNA" \
            --fastx \
            --threads "${THREADS}"
        echo "  ✓ rRNA removed"
    fi

    # ── Step 3: STAR alignment ────────────────────────────────────────────────
    RAW_BAM="${PHASE3}/03_aligned/${sample_name}_Aligned.sortedByCoord.out.bam"
    if [ ! -f "${RAW_BAM}" ] && [ ! -f "${RPF_BAM}" ]; then
        echo "  STAR alignment..."
        STAR \
            --runThreadN "${THREADS}" \
            --genomeDir "${STAR_INDEX}" \
            --readFilesIn "${NORRNA}" \
            --readFilesCommand zcat \
            --sjdbGTFfile "${GTF}" \
            --outFileNamePrefix "${PHASE3}/03_aligned/${sample_name}_" \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes NH HI AS NM MD \
            --outFilterMultimapNmax 1 \
            --outFilterMismatchNoverReadLmax 0.04 \
            --outFilterMatchNmin 20 \
            --limitBAMsortRAM 20000000000
        samtools index "${RAW_BAM}"
        echo "  ✓ Aligned"
    fi

    # ── Step 4: RPF length filter ─────────────────────────────────────────────
    if [ ! -f "${RPF_BAM}" ]; then
        echo "  Filtering RPF lengths (${RPF_MIN}–${RPF_MAX} nt)..."
        samtools view -h "${RAW_BAM}" | \
        awk -v min="${RPF_MIN}" -v max="${RPF_MAX}" '
        BEGIN{OFS="\t"}
        /^@/ {print; next}
        { if (length($10) >= min && length($10) <= max) print }' | \
        samtools sort -@ "${THREADS}" -o "${RPF_BAM}"
        samtools index "${RPF_BAM}"
        echo "  ✓ RPF BAM written"
    fi

    # ── Step 5: P-site assignment ─────────────────────────────────────────────
    if [ ! -f "${PSITE_BED}" ]; then
        echo "  P-site assignment..."
        samtools view "${RPF_BAM}" | \
        awk 'BEGIN{OFS="\t"}
        {
            chrom   = $3
            start   = $4 - 1
            strand  = and($2, 16) ? "-" : "+"
            readlen = length($10)
            if      (readlen == 28 || readlen == 29) offset = 12
            else if (readlen == 30 || readlen == 31) offset = 13
            else                                     offset = 12
            if (strand == "+")
                psite = start + offset
            else
                psite = start + readlen - offset - 1
            print chrom, psite, psite+1, ".", ".", strand
        }' > "${PSITE_BED}"
        echo "  ✓ P-sites written: ${PSITE_BED}"
    fi

done < "${SAMPLES}"

echo ""
echo "============================================================"
echo "  PHASE 3 COMPLETE — $(date)"
echo "============================================================"
