# Test Dataset — HeLa Asynchronous (Validation Run)

Two validation tests are provided depending on available time and resources.

---

## Quick Test (~10 minutes) — Recommended

Uses precomputed alignment files from Zenodo. Runs Phases 4+5 only.

**Zenodo DOI: https://doi.org/10.5281/zenodo.20084600**

### Prerequisites
```bash
# 1. Install environment
conda env create -f environment.yml
conda activate ribowin

# 2. Download GTF
wget -P data/raw/ https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz
gunzip data/raw/Homo_sapiens.GRCh38.115.gtf.gz
```

### Run
```bash
bash test/run_test_quick.sh
```

Zenodo data downloads automatically. No manual file preparation needed.

### Expected output
- **12 concordant sORFs** in `results/phase4/ribo_HeLa_async_common_translated_orfs.csv`
- **RPL26P19** present — pseudogene-derived microprotein candidate

### What Zenodo provides
| File | Description |
|------|-------------|
| `stage1_cleaned_sorfs.csv` | sORF table from GRCh38 transcriptome scan |
| `HeLa_async_RNA_rep1/2_sorted.bam` | STAR-aligned RNA-seq BAMs + indices |
| `HeLa_async_RIBO_rep1/2_psites.bed` | P-site BEDs with empirical offsets |

---

## Full Pipeline Test (~5-6 hours)

Downloads raw FASTQ data from SRA and runs Phases 2-5 from scratch.
Validates end-to-end reproducibility including alignment and P-site calibration.

### Prerequisites
```bash
# 1. Install environment
conda env create -f environment.yml
conda activate ribowin

# 2. All reference files in data/raw/ (see README.md)
#    Genome FASTA, GTF, rRNA databases

# 3. Precomputed Phase 1 output from GitHub Releases
mkdir -p results/phase1
# Download stage1_cleaned_sorfs.csv from Releases page
# Place at: results/phase1/stage1_cleaned_sorfs.csv
```

### Run
```bash
bash test/run_test_full.sh
```

Raw data is downloaded automatically from SRA:
- SRR3306581, SRR3306582 (HeLa async RNA-seq)
- SRR3306588, SRR3306589 (HeLa async Ribo-seq)

### Expected output
- **12 concordant sORFs** reproducible across both replicates
- **RPL26P19** present — consistent with quick test and paper results

---

## Dataset

HeLa asynchronous cells, 2 biological replicates.
Source: Aviner et al. (2017), *Nature Structural & Molecular Biology*, GSE79664.
