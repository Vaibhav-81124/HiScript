#!/usr/bin/env bash
# =============================================================================
#  run_phase3.sh — Ribo-seq processing
#
#  Runs for each ribo sample in samples.tsv:
#    1. SRA download (if SRR provided)
#    2. Trimmomatic trim + size-select (RPF_MIN-RPF_MAX)
#    3. rRNA removal (SortMeRNA)
#    4. STAR alignment
#    5. RPF length filter (samtools)
#    6a. Empirical P-site calibration
#    6b. P-site assignment -> BED
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
LEADING=$(read_cfg riboseq.trimmomatic_leading)
TRAILING=$(read_cfg riboseq.trimmomatic_trailing)
OUTDIR=$(read_cfg output_dir)

RIBO_DATA="data/ribo_seq"
mkdir -p "${RIBO_DATA}"

PHASE3="${OUTDIR}/phase3"
mkdir -p "${PHASE3}/01_trimmed" \
         "${PHASE3}/02_norrna" \
         "${PHASE3}/03_aligned" \
         "${PHASE3}/04_psites"

echo "============================================================"
echo "  PHASE 3 -- Ribo-seq Pipeline (Trimmomatic)"
echo "  Config  : ${CONFIG}"
echo "  Samples : ${SAMPLES}"
echo "  $(date)"
echo "============================================================"

if [ ! -f "${STAR_INDEX}/SA" ]; then
    echo "ERROR: STAR index not found at ${STAR_INDEX}. Run phase 2 first."
    exit 1
fi

while IFS=$'\t' read -r sample_name cell_type data_type layout replicate srr fastq_r1 fastq_r2; do
    [[ "${sample_name}" =~ ^#.*$ || "${sample_name}" == "sample_name" ]] && continue
    [[ "${data_type}" != "ribo" ]] && continue

    echo ""
    echo "------------------------------------------------------------"
    echo "  Sample: ${sample_name}"
    echo "------------------------------------------------------------"

    # -- Download if SRR provided ---------------------------------------------
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

    # -- Step 1: Trimmomatic + size select ------------------------------------
    # Ribo-seq: trim adapter, then size-select RPF_MIN-RPF_MAX nt in one pass.
    # MINLEN/MAXLEN enforced via MINLEN and a post-trim awk filter since
    # Trimmomatic does not have MAXLEN; we filter with awk after trimming.
    if [ ! -f "${TRIMMED}" ]; then
        echo "  Trimmomatic (single-end, Ribo-seq)..."

        ADAPTER_FA="${PHASE3}/01_trimmed/${sample_name}_adapter.fa"
        printf ">adapter\n%s\n" "${ADAPTER}" > "${ADAPTER_FA}"

        TRIMMED_FULL="${PHASE3}/01_trimmed/${sample_name}_trimmed_full.fastq.gz"

        trimmomatic SE \
            -threads "${THREADS}" \
            -phred33 \
            "${fastq_r1}" \
            "${TRIMMED_FULL}" \
            "ILLUMINACLIP:${ADAPTER_FA}:2:30:10" \
            "LEADING:${LEADING}" \
            "TRAILING:${TRAILING}" \
            "MINLEN:${RPF_MIN}" \
            2>&1 | tee "${PHASE3}/01_trimmed/${sample_name}_trimmomatic.log"

        # Size-select: keep only RPF_MIN to RPF_MAX nt reads
        echo "  Size-selecting ${RPF_MIN}-${RPF_MAX} nt reads..."
        zcat "${TRIMMED_FULL}" | \
        awk -v min="${RPF_MIN}" -v max="${RPF_MAX}" '
        NR%4==1 { header=$0 }
        NR%4==2 { seq=$0 }
        NR%4==3 { plus=$0 }
        NR%4==0 {
            if (length(seq) >= min && length(seq) <= max)
                print header"\n"seq"\n"plus"\n"$0
        }' | gzip > "${TRIMMED}"

        rm -f "${TRIMMED_FULL}"
        echo "  Trimmed and size-selected"
    fi

    # -- Step 2: rRNA removal -------------------------------------------------
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
        echo "  rRNA removed"
    fi

    # -- Step 3: STAR alignment -----------------------------------------------
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
            --outFilterMismatchNmax 2 \
            --outFilterMatchNmin 20 \
            --alignIntronMax 1000000 \
            --alignMatesGapMax 1000000 \
            --limitBAMsortRAM 20000000000
        samtools index "${RAW_BAM}"
        echo "  Aligned"
    fi

    # -- Step 4: RPF length filter --------------------------------------------
    if [ ! -f "${RPF_BAM}" ]; then
        echo "  Filtering RPF lengths (${RPF_MIN}-${RPF_MAX} nt)..."
        samtools view -h "${RAW_BAM}" | \
        awk -v min="${RPF_MIN}" -v max="${RPF_MAX}" '
        BEGIN{OFS="\t"}
        /^@/ {print; next}
        { if (length($10) >= min && length($10) <= max) print }' | \
        samtools sort -@ "${THREADS}" -o "${RPF_BAM}"
        samtools index "${RPF_BAM}"
        echo "  RPF BAM written"
    fi

<<<<<<< HEAD
    # -- Step 5a: Empirical P-site calibration --------------------------------
=======
    # ── Step 5: P-site assignment ─────────────────────────────────────────────
    # Step 5a: Empirical P-site calibration
>>>>>>> b38740268b6318bad29af99828006001b49a86d7
    OFFSETS_JSON="${PHASE3}/04_psites/${sample_name}_psite_offsets.json"

    if [ ! -f "${OFFSETS_JSON}" ]; then
        echo "  Calibrating P-site offsets empirically..."
        python3 scripts/psite_calibration.py \
            --bam      "${RPF_BAM}" \
            --gtf      "${GTF}" \
            --sample   "${sample_name}" \
            --rpf_min  "${RPF_MIN}" \
            --rpf_max  "${RPF_MAX}" \
            --outdir   "${PHASE3}/04_psites"
<<<<<<< HEAD
        echo "  Offsets calibrated"
    fi

    # -- Step 5b: P-site assignment -------------------------------------------
    if [ ! -f "${PSITE_BED}" ]; then
        echo "  Assigning P-sites..."
        python3 scripts/psite_assignment.py \
            --bam      "${RPF_BAM}" \
            --offsets  "${OFFSETS_JSON}" \
            --sample   "${sample_name}" \
            --rpf_min  "${RPF_MIN}" \
            --rpf_max  "${RPF_MAX}" \
            --outdir   "${PHASE3}/04_psites"
        echo "  P-sites written: ${PSITE_BED}"
=======
        echo "  Offsets calibrated: ${OFFSETS_JSON}"
>>>>>>> b38740268b6318bad29af99828006001b49a86d7
    fi

    # Step 5b: P-site assignment using calibrated offsets
    if [ ! -f "${PSITE_BED}" ]; then
        echo "  Assigning P-sites..."
        python3 scripts/psite_assignment.py \
            --bam      "${RPF_BAM}" \
            --offsets  "${OFFSETS_JSON}" \
            --sample   "${sample_name}" \
            --rpf_min  "${RPF_MIN}" \
            --rpf_max  "${RPF_MAX}" \
            --outdir   "${PHASE3}/04_psites"
        echo "  P-sites written: ${PSITE_BED}"
    fi
done < "${SAMPLES}"

echo ""
echo "============================================================"
echo "  PHASE 3 COMPLETE -- $(date)"
echo "============================================================"
