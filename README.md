# HiScript Tool

**Computational discovery and translational validation of small open reading frames (sORFs) and microproteins from multi-omics data.**

---

## Overview

This pipeline integrates transcriptomic sORF scanning with ribosome profiling and RNA-seq to identify and validate translated microproteins. It was developed using HeLa and other human cell-line data aligned to GRCh38.

```
Phase 1 → sORF Discovery          (Python)
Phase 2 → RNA-seq Alignment        (Bash / STAR)
Phase 3 → Ribo-seq Processing      (Bash / STAR / SortMeRNA)
Phase 4 → Translation Evidence     (Python + bedtools)
Phase 5 → Translation Efficiency   (Python)
```

---

## Quick Start

### 1. Install dependencies

```bash
conda env create -f environment.yml
conda activate sorf-tool
```

### 2. Place reference files

Download the following into `data/raw/`:

| File | Source |
|------|--------|
| `Homo_sapiens.GRCh38.dna.primary_assembly.fa` | [Ensembl](https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/) |
| `Homo_sapiens.GRCh38.cdna.all.fa` | [Ensembl](https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/cdna/) |
| `Homo_sapiens.GRCh38.115.gtf` | [Ensembl](https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/) |
| `rRNA_db/` (SortMeRNA databases) | [SortMeRNA GitHub](https://github.com/sortmerna/sortmerna/tree/master/data/rRNA_databases) |

Update paths in `config/config.yaml` if you place files elsewhere.

### 3. Configure your samples

Edit `config/samples.tsv` to list your samples:

```
sample_name         cell_type  data_type  layout  replicate  srr         fastq_r1  fastq_r2
HeLa_M_RNA_rep1     HeLa_M     rna        single  rep1       SRR3306577
HeLa_M_RIBO_rep1    HeLa_M     ribo       single  rep1       SRR3306585
```

- Fill `srr` to auto-download from SRA, **or** fill `fastq_r1`/`fastq_r2` for local files.
- For paired-end RNA-seq, set `layout=paired` and provide `fastq_r2`.

### 4. Run

**Full pipeline:**
```bash
bash run_all.sh
```

**Skip sORF discovery (use precomputed `stage1_novel_sorfs.csv`):**
```bash
bash run_all.sh --skip_discovery
```

**Run specific phases only:**
```bash
bash run_all.sh --start_phase 2 --end_phase 3
```

**Run a single phase:**
```bash
bash phase1_discovery/run_phase1.sh
bash phase2_rnaseq/run_phase2.sh
bash phase3_riboseq/run_phase3.sh
bash phase4_translation/run_phase4.sh
bash phase5_TE/run_phase5.sh
```

---

## Repository Structure

```
sorf-tool/
├── config/
│   ├── config.yaml             # All parameters — edit this
│   └── samples.tsv             # Sample manifest — edit this
├── data/
│   ├── raw/                    # Reference files (user-supplied)
│   └── processed/              # Auto-created intermediate files
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
│   ├── sorf_discovery.py       # Phase 1a: scan transcriptome
│   ├── sorf_cleaning.py        # Phase 1b: dedup + filter
│   ├── make_bed.py             # sORF CSV → genomic BED
│   ├── make_cds_bed.py         # GTF → CDS BED
│   ├── stage4_translation.py   # Ribo count merge + CDS filter
│   ├── stage5_periodicity.py   # Triplet periodicity filter
│   ├── merge_reps.py           # Intersect translated ORFs across reps
│   ├── TE_compute.py           # Compute Translation Efficiency
│   └── TE_filter.py            # Filter TE to translated ORFs
├── results/                    # Auto-created outputs
├── run_all.sh                  # Master runner
├── environment.yml             # Conda environment
└── README.md
```

---

## Output Files

| Phase | Key outputs |
|-------|-------------|
| Phase 1 | `results/phase1/stage1_cleaned_sorfs.csv`, `sorfs_genomic.bed`, `cds.bed` |
| Phase 2 | `results/phase2/03_aligned/<sample>_sorted.bam` |
| Phase 3 | `results/phase3/04_psites/<sample>_psites.bed` |
| Phase 4 | `results/phase4/ribo_<cell>_common_translated_orfs.csv` |
| Phase 5 | `results/phase5/<cell>_translated_orfs_filtered_withTE.csv` |

---

## Parameters

All parameters are in `config/config.yaml`. Key ones:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `discovery.min_aa` | 10 | Min peptide length for scanning |
| `discovery.min_aa_after_cleaning` | 15 | Min length after deduplication |
| `riboseq.rpf_min` / `rpf_max` | 26 / 34 | RPF size selection window (nt) |
| `translation.min_ribo_reads` | 10 | Min Ribo reads to call translated |
| `translation.frame0_threshold` | 0.55 | Min frame-0 fraction for periodicity |
| `te.pseudocount` | 1 | Added to RNA denominator to avoid div/0 |

---

## Notes on Phase 1

sORF discovery (`sorf_discovery.py`) requires the full GRCh38 cDNA FASTA and GTF (~5 GB) and takes several hours. A precomputed `stage1_novel_sorfs.csv` is available in [Releases](../../releases) — use `--skip_discovery` to use it and skip to cleaning.

---

## Citation

> *[Your paper title here]* — [Authors], [Journal], [Year]

---

## Contact

[Your name / lab / GitHub issues link]
