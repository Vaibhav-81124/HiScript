# RiboWin

**Computational discovery and translational validation of small open reading frames (sORFs) and microproteins from multi-omics Ribo-seq and RNA-seq data.**

---

## Overview

RiboWin is a modular, reproducible pipeline that integrates transcriptome-wide sORF scanning with ribosome profiling (Ribo-seq) and RNA-seq to identify and validate translated microproteins. It was developed and benchmarked on human cell lines (HeLa, keratinocyte, fibroblast) aligned to GRCh38.

```
Phase 1 → sORF Discovery           scan transcriptome for candidate ORFs
Phase 2 → RNA-seq Alignment         Cutadapt + STAR
Phase 3 → Ribo-seq Processing       Cutadapt + SortMeRNA + STAR + P-site calibration
Phase 4 → Translation Evidence      Ribo counts + periodicity filter + replicate concordance
Phase 5 → Translation Efficiency    TE calculation per cell type
```

Key features:
- Zero hardcoded paths — all parameters in `config/config.yaml`
- Auto-downloads FASTQs from SRA given accession numbers
- Empirical P-site offset calibration per sample (no hardcoded offsets)
- Mandatory two-replicate concordance filter
- Handles single-end and paired-end RNA-seq

---

## Requirements

- Linux or macOS (WSL2 on Windows)
- conda / mamba
- ~50 GB disk space for reference files + results
- 16+ GB RAM recommended for STAR alignment

---

## Installation

```bash
git clone https://github.com/Vaibhav-81124/RiboWin.git
cd RiboWin
conda env create -f environment.yml
conda activate ribowin
```

---

## Reference Files

Download the following into `data/raw/`. Exact filenames matter — they match the defaults in `config/config.yaml`.

**Genome FASTA:**
```bash
wget -P data/raw/ https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip data/raw/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
```

**GTF:**
```bash
wget -P data/raw/ https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz
gunzip data/raw/Homo_sapiens.GRCh38.115.gtf.gz
```

**SortMeRNA rRNA databases:**
```bash
mkdir -p data/raw/rRNA_db
wget -P data/raw/rRNA_db/ https://github.com/sortmerna/sortmerna/raw/master/data/rRNA_databases/silva-euk-28s-id98.fasta
wget -P data/raw/rRNA_db/ https://github.com/sortmerna/sortmerna/raw/master/data/rRNA_databases/silva-euk-18s-id95.fasta
wget -P data/raw/rRNA_db/ https://github.com/sortmerna/sortmerna/raw/master/data/rRNA_databases/rfam-5s-database-id98.fasta
wget -P data/raw/rRNA_db/ https://github.com/sortmerna/sortmerna/raw/master/data/rRNA_databases/rfam-5.8s-database-id98.fasta
```

---

## Quick Start — Validation Test

Two tests are provided. Run the quick test first to verify your setup.

### Quick test (~10 minutes)
Uses precomputed alignment files from Zenodo (DOI: https://doi.org/10.5281/zenodo.20084600).
Runs Phases 4+5 only. Zenodo data downloads automatically.

```bash
# Only GTF needed — no genome FASTA required for quick test
bash test/run_test_quick.sh
```
or 
If you wish to run phase 1, download cDNA from the link below
[Releases](https://github.com/Vaibhav-81124/RiboWin/releases)

Expected: **12 concordant sORFs including RPL26P19**

### Full pipeline test (~5-6 hours)
Downloads raw FASTQs from SRA and runs Phases 2-5 from scratch.
Requires all reference files and precomputed Phase 1 CSV from Releases.

```bash
# Download stage1_cleaned_sorfs.csv from Releases -> results/phase1/
bash test/run_test_full.sh
```

See `test/README_test.md` for full details on both tests.

---

## Running on Your Own Data

### Step 1 — Configure samples

Edit `config/samples.tsv` with your sample information:

```
sample_name     cell_type   data_type  layout  replicate  srr
MY_RNA_rep1     MY_CELL     rna        single  rep1       SRR0000001
MY_RNA_rep2     MY_CELL     rna        single  rep2       SRR0000002
MY_RIBO_rep1    MY_CELL     ribo       single  rep1       SRR0000003
MY_RIBO_rep2    MY_CELL     ribo       single  rep2       SRR0000004
```

- `srr`: SRA accession — pipeline downloads automatically via `prefetch` + `fasterq-dump`
- `fastq_r1` / `fastq_r2`: use these instead of `srr` if you have local FASTQ files
- `layout`: `single` or `paired` — affects Cutadapt and featureCounts settings
- Each `cell_type` must have **at least 2 replicates** for both `rna` and `ribo`

### Step 2 — Edit parameters (optional)

All parameters are in `config/config.yaml`. Defaults work for human GRCh38 data.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `threads` | 12 | CPU threads for all tools |
| `discovery.min_aa` | 10 | Min ORF length for scanning (aa) |
| `discovery.min_aa_after_cleaning` | 15 | Min length after deduplication |
| `riboseq.rpf_min` / `rpf_max` | 26 / 34 | RPF size selection window (nt) |
| `translation.min_ribo_reads` | 10 | Min Ribo-seq reads to call translated |
| `translation.frame0_threshold` | 0.55 | Min frame-0 fraction for periodicity |
| `te.pseudocount` | 1 | Pseudocount added to RNA denominator |

### Step 3 — Get precomputed Phase 1 output

Phase 1 (sORF scanning) requires the full cDNA FASTA.
Download our precomputed output instead:[Releases](https://github.com/Vaibhav-81124/RiboWin/releases)

```bash
mkdir -p results/phase1
# Download stage1_cleaned_sorfs.csv from GitHub Releases
# Place at: results/phase1/stage1_cleaned_sorfs.csv
```

Then run:
```bash
bash run_all.sh --skip_discovery
```

To run Phase 1 from scratch (only needed for non-human genomes or custom parameters):
```bash
bash run_all.sh
```

---

## Running the Pipeline

**Full pipeline:**
```bash
bash run_all.sh
```

**Skip Phase 1 (recommended):**
```bash
bash run_all.sh --skip_discovery
```

**Specific phase range:**
```bash
bash run_all.sh --start_phase 2 --end_phase 3
```

**Individual phases:**
```bash
bash phase1_discovery/run_phase1.sh
bash phase2_rnaseq/run_phase2.sh
bash phase3_riboseq/run_phase3.sh
bash phase4_translation/run_phase4.sh
bash phase5_TE/run_phase5.sh
```

All phases are idempotent — if an output file already exists the step is skipped. Safe to re-run after a crash.

---

## Full Reproduction from Raw Data

The complete pipeline can be reproduced from raw SRA accessions:

```bash
# Edit config/samples.tsv with your SRR accessions
bash run_all.sh --skip_discovery
```

**Tested on:** Ubuntu 22.04, 32 cores, 64 GB RAM
**Runtime:** ~5-6 hours for 2 replicates (RNA-seq + Ribo-seq)

Raw data for the HeLa asynchronous condition used in the paper:

| Sample | SRR | Type |
|--------|-----|------|
| HeLa async RNA rep1 | SRR3306581 | RNA-seq |
| HeLa async RNA rep2 | SRR3306582 | RNA-seq |
| HeLa async Ribo rep1 | SRR3306588 | Ribo-seq |
| HeLa async Ribo rep2 | SRR3306589 | Ribo-seq |

GEO accession: GSE79664 (Aviner et al. 2017)

---

## Output Files

| Phase | Key output | Description |
|-------|------------|-------------|
| Phase 1 | `results/phase1/stage1_cleaned_sorfs.csv` | Deduplicated sORF table |
| Phase 1 | `results/phase1/sorfs_genomic.bed` | Genomic coordinates BED |
| Phase 1 | `results/phase1/cds.bed` | CDS regions from GTF |
| Phase 2 | `results/phase2/03_aligned/<sample>_sorted.bam` | Aligned RNA-seq BAM |
| Phase 3 | `results/phase3/04_psites/<sample>_psites.bed` | P-site positions BED |
| Phase 3 | `results/phase3/04_psites/<sample>_psite_offsets.json` | Empirical P-site offsets |
| Phase 4 | `results/phase4/ribo_<cell>_common_translated_orfs.csv` | Concordant translated ORFs |
| Phase 5 | `results/phase5/<cell>_translated_orfs_filtered_withTE.csv` | Final ORFs with TE scores |

---

## Repository Structure

```
RiboWin/
├── config/
│   ├── config.yaml             # All parameters — edit this
│   └── samples.tsv             # Sample manifest — edit this
├── data/
│   └── raw/                    # Place reference files here
├── phase1_discovery/
│   └── run_phase1.sh
├── phase2_rnaseq/
│   └── run_phase2.sh
├── phase3_riboseq/
│   └── run_phase3.sh
├── phase4_translation/
│   └── run_phase4.sh
├── phase5_TE/
│   └── run_phase5.sh
├── scripts/
│   ├── sorf_discovery.py       # Phase 1a: scan transcriptome for sORFs
│   ├── sorf_cleaning.py        # Phase 1b: deduplicate and filter
│   ├── make_bed.py             # sORF CSV to genomic BED
│   ├── make_cds_bed.py         # GTF to CDS BED
│   ├── psite_calibration.py    # Empirical P-site offset calibration
│   ├── psite_assignment.py     # P-site assignment using calibrated offsets
│   ├── stage4_translation.py   # Ribo count merge + CDS overlap removal
│   ├── stage5_periodicity.py   # Triplet periodicity filter
│   ├── merge_reps.py           # Concordance filter across replicates
│   ├── TE_compute.py           # Translation efficiency calculation
│   └── TE_filter.py            # Filter TE table to translated ORFs
├── test/
│   ├── samples_test.tsv        # HeLa async test sample manifest
│   ├── run_test_quick.sh       # Quick test: Zenodo data, Phases 4+5, ~10 min
│   ├── run_test_full.sh        # Full test: SRA download, Phases 2-5, ~5-6 hrs
│   └── README_test.md          # Test instructions and expected output
├── run_all.sh                  # Master pipeline runner
├── environment.yml             # Conda environment
└── README.md
```

---

## Data Availability

- **Code**: https://github.com/Vaibhav-81124/RiboWin (MIT License)
- **Precomputed test data**: https://doi.org/10.5281/zenodo.20084600

