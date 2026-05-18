#!/usr/bin/env bash
# =============================================================================
#  run_phase3.sh -- Ribo-seq processing 
#
#  For each ribo sample in samples.tsv:
#    1. SRA download (if SRR provided)
#    2. Cutadapt trim + size-select
#    3. rRNA removal (SortMeRNA)
#    4. STAR alignment
#    5. RPF length filtering
#    6. Fixed-offset P-site assignment
#
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

# =============================================================================
# CONFIG READER
# =============================================================================

read_cfg() {
python3 - "$CONFIG" "$1" << 'PYEOF'
import yaml
import sys

config = sys.argv[1]
key = sys.argv[2]

with open(config) as f:
    c = yaml.safe_load(f)

v = c
for k in key.split('.'):
    v = v[k]

print(v)
PYEOF
}

# =============================================================================
# LOAD CONFIG
# =============================================================================

GTF=$(read_cfg reference.gtf)
STAR_INDEX=$(read_cfg reference.star_index)
RRNA_DB=$(read_cfg reference.rrna_db)

THREADS=$(read_cfg threads)

ADAPTER=$(read_cfg riboseq.adapter)
RPF_MIN=$(read_cfg riboseq.rpf_min)
RPF_MAX=$(read_cfg riboseq.rpf_max)

OUTDIR=$(read_cfg output_dir)

RIBO_DATA="data/ribo_seq"
mkdir -p "${RIBO_DATA}"

PHASE3="${OUTDIR}/phase3"

mkdir -p \
    "${PHASE3}/01_trimmed" \
    "${PHASE3}/02_norrna" \
    "${PHASE3}/03_aligned" \
    "${PHASE3}/04_psites"

echo "============================================================"
echo "  PHASE 3 -- Ribo-seq Pipeline"
echo "  Manuscript-consistent fixed-offset version"
echo "  Config  : ${CONFIG}"
echo "  Samples : ${SAMPLES}"
echo "  $(date)"
echo "============================================================"

# =============================================================================
# SANITY CHECK
# =============================================================================

if [ ! -f "${STAR_INDEX}/SA" ]; then
    echo "ERROR: STAR index not found at ${STAR_INDEX}"
    exit 1
fi

# =============================================================================
# PROCESS SAMPLES
# =============================================================================

while IFS=$'\t' read -r sample_name cell_type data_type layout replicate srr fastq_r1 fastq_r2 || [[ -n "$sample_name" ]]; do

    sample_name=$(echo "$sample_name" | tr -d '\r')
    data_type=$(echo "$data_type" | tr -d '\r')
    srr=$(echo "$srr" | tr -d '\r')

    [[ "${sample_name}" =~ ^#.*$ || "${sample_name}" == "sample_name" ]] && continue
    [[ "${data_type}" != "ribo" ]] && continue

    echo ""
    echo "------------------------------------------------------------"
    echo "  Sample: ${sample_name}"
    echo "------------------------------------------------------------"

    # =========================================================================
    # DOWNLOAD
    # =========================================================================

    if [ -n "${srr}" ]; then

        FQ="${RIBO_DATA}/${srr}.fastq.gz"

        if [ ! -f "${FQ}" ]; then

            echo "  Downloading ${srr}..."

            prefetch "${srr}" --output-directory "${RIBO_DATA}"

            fasterq-dump \
                "${RIBO_DATA}/${srr}" \
                --outdir "${RIBO_DATA}" \
                --threads "${THREADS}"

            if [ -f "${RIBO_DATA}/${srr}.fastq" ]; then

                gzip "${RIBO_DATA}/${srr}.fastq"
                fastq_r1="${RIBO_DATA}/${srr}.fastq.gz"

            elif [ -f "${RIBO_DATA}/${srr}_1.fastq" ]; then

                gzip "${RIBO_DATA}/${srr}_1.fastq"
                fastq_r1="${RIBO_DATA}/${srr}_1.fastq.gz"

            else
                echo "ERROR: FASTQ not found for ${srr}"
                continue
            fi

        else
            fastq_r1="${FQ}"
        fi
    fi

    # =========================================================================
    # PATHS
    # =========================================================================

    TRIMMED="${PHASE3}/01_trimmed/${sample_name}_trimmed.fastq.gz"

    NORRNA_PREFIX="${PHASE3}/02_norrna/${sample_name}_norRNA"
    NORRNA="${NORRNA_PREFIX}.fq.gz"

    STAR_PREFIX="${PHASE3}/03_aligned/${sample_name}_"

    RAW_BAM="${STAR_PREFIX}Aligned.sortedByCoord.out.bam"

    RPF_BAM="${PHASE3}/03_aligned/${sample_name}_RPFs.bam"

    PSITE_BED="${PHASE3}/04_psites/${sample_name}_psites.bed"

    # =========================================================================
    # STEP 1: CUTADAPT
    # =========================================================================

    if [ ! -f "${TRIMMED}" ]; then

        echo "  Cutadapt trimming + size selection..."

        cutadapt \
            -a "${ADAPTER}" \
            --quality-cutoff 20 \
            --minimum-length "${RPF_MIN}" \
            --maximum-length "${RPF_MAX}" \
            --cores "${THREADS}" \
            -o "${TRIMMED}" \
            "${fastq_r1}" \
            2>&1 | tee "${PHASE3}/01_trimmed/${sample_name}_cutadapt.log"

        echo "  Trimmed and size-selected"
    fi

    # =========================================================================
    # STEP 2: rRNA REMOVAL
    # =========================================================================

    if [ ! -f "${NORRNA}" ]; then

        echo "  rRNA removal (SortMeRNA)..."

        WORKDIR="${PHASE3}/02_norrna/sortmerna_${sample_name}_workdir"

        rm -rf "${WORKDIR}"
        mkdir -p "${WORKDIR}"

        sortmerna \
            --ref "${RRNA_DB}/silva-euk-28s-id98.fasta" \
            --ref "${RRNA_DB}/silva-euk-18s-id95.fasta" \
            --ref "${RRNA_DB}/rfam-5s-database-id98.fasta" \
            --ref "${RRNA_DB}/rfam-5.8s-database-id98.fasta" \
            --reads "${TRIMMED}" \
            --workdir "${WORKDIR}" \
            --other "${NORRNA_PREFIX}" \
            --fastx \
            --threads "${THREADS}"

        echo "  rRNA removed"
    fi

    # =========================================================================
    # STEP 3: STAR ALIGNMENT
    # =========================================================================

    if [ ! -f "${RAW_BAM}" ]; then

        echo "  STAR alignment..."

        STAR \
            --runThreadN "${THREADS}" \
            --genomeDir "${STAR_INDEX}" \
            --readFilesIn "${NORRNA}" \
            --readFilesCommand zcat \
            --sjdbGTFfile "${GTF}" \
            --outFileNamePrefix "${STAR_PREFIX}" \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes NH HI AS NM MD \
            --outFilterMultimapNmax 1 \
            --outFilterMismatchNoverReadLmax 0.04 \
            --outFilterMatchNmin 20 \
            --limitBAMsortRAM 20000000000 \
            --quantMode TranscriptomeSAM

        samtools index "${RAW_BAM}"

        echo "  Alignment complete"
    fi

    # =========================================================================
    # STEP 4: FILTER RPF LENGTHS
    # =========================================================================

    if [ ! -f "${RPF_BAM}" ]; then

        echo "  Filtering RPF lengths (${RPF_MIN}-${RPF_MAX} nt)..."

        samtools view -h "${RAW_BAM}" | \
        awk -v min="${RPF_MIN}" -v max="${RPF_MAX}" '
        BEGIN{OFS="\t"}
        /^@/ {print; next}
        {
            if (length($10) >= min && length($10) <= max)
                print
        }' | \
        samtools sort -@ "${THREADS}" -o "${RPF_BAM}"

        samtools index "${RPF_BAM}"

        echo "  RPF BAM written"
    fi

    # =========================================================================
    # STEP 5: P-SITE ASSIGNMENT
    # =========================================================================

    if [ ! -f "${PSITE_BED}" ]; then

        echo "  Assigning fixed-offset P-sites..."

        samtools view "${RPF_BAM}" | \
        awk 'BEGIN{OFS="\t"}
        {
            chrom=$3
            start=$4-1
            strand = and($2,16) ? "-" : "+"
            readlen=length($10)

            if (readlen==28 || readlen==29)
                offset=12
            else if (readlen==30 || readlen==31)
                offset=13
            else
                offset=12

            if (strand=="+")
                psite=start+offset
            else
                psite=start+readlen-offset-1

            print chrom, psite, psite+1, ".", ".", strand
        }' > "${PSITE_BED}"

        echo "  P-sites written: ${PSITE_BED}"
    fi

done < "${SAMPLES}"

echo ""
echo "============================================================"
echo "  PHASE 3 COMPLETE -- $(date)"
echo "============================================================"
