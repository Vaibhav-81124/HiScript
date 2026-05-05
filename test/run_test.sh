#!/usr/bin/env bash
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

# 🔴 Fix CRLF issues globally (safe no-op if already clean)
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
echo "  PHASE 3 -- Ribo-seq Pipeline"
echo "============================================================"

if [ ! -f "${STAR_INDEX}/SA" ]; then
    echo "ERROR: STAR index not found"
    exit 1
fi

# ✅ FIXED LOOP (handles last line + CRLF safely)
while IFS=$'\t' read -r sample_name cell_type data_type layout replicate srr fastq_r1 fastq_r2 || [[ -n "$sample_name" ]]; do

    # 🔴 sanitize fields (CRLF-safe)
    sample_name=$(echo "$sample_name" | tr -d '\r')
    data_type=$(echo "$data_type" | tr -d '\r')
    srr=$(echo "$srr" | tr -d '\r')

    [[ "${sample_name}" =~ ^#.*$ || "${sample_name}" == "sample_name" ]] && continue
    [[ "${data_type}" != "ribo" ]] && continue

    echo ""
    echo "------------------------------------------------------------"
    echo "  Sample: ${sample_name} (${srr})"
    echo "------------------------------------------------------------"

    # ✅ Robust download block
    if [ -n "${srr}" ]; then
        FQ="${RIBO_DATA}/${srr}.fastq.gz"

        if [ ! -f "${FQ}" ]; then
            echo "  Downloading ${srr}..."

            prefetch "${srr}" --output-directory "${RIBO_DATA}"

            fasterq-dump "${RIBO_DATA}/${srr}" \
                --outdir "${RIBO_DATA}" \
                --threads "${THREADS}"

            # Handle possible outputs
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

    # ---- rest of your pipeline unchanged ----

done < "${SAMPLES}"

echo "============================================================"
echo "  PHASE 3 COMPLETE"
echo "============================================================"