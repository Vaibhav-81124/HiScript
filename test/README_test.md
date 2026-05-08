# Test Dataset — HeLa Asynchronous (Validation Run)

Validates the pipeline using precomputed HeLa asynchronous data hosted on Zenodo.
Only Phases 4 and 5 are run (translation evidence + TE). No alignment needed.

**Runtime: ~10 minutes**

---

## Dataset

HeLa asynchronous cells, 2 biological replicates.
Source: Aviner et al. (2017), *Nature Structural & Molecular Biology*, GSE79664.

Precomputed files hosted on Zenodo: https://doi.org/10.5281/zenodo.20084600

---

## Setup

**Step 1 — Install environment:**
```bash
conda env create -f environment.yml
conda activate ribowin
```

**Step 2 — Download GTF into `data/raw/`:**
```bash
wget -P data/raw/ https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz
gunzip data/raw/Homo_sapiens.GRCh38.115.gtf.gz
```

**Step 3 — Download precomputed files from Zenodo:**

Extract the Zenodo archive into the repo root. It will populate:
```
results/
├── phase1/
│   └── stage1_cleaned_sorfs.csv
├── phase2/
│   └── 03_aligned/
│       ├── HeLa_async_RNA_rep1_sorted.bam  + .bai
│       └── HeLa_async_RNA_rep2_sorted.bam  + .bai
└── phase3/
    └── 04_psites/
        ├── HeLa_async_RIBO_rep1_psites.bed
        └── HeLa_async_RIBO_rep2_psites.bed
```

**Step 4 — Run:**
```bash
bash test/run_test.sh
```

---

## Expected output

`results/phase4/ribo_HeLa_async_common_translated_orfs.csv`

- **12 concordant sORFs** reproducible across both replicates
- **RPL26P19** present on chr5 — pseudogene-derived microprotein candidate

The test script checks both conditions and prints PASS/FAIL automatically.

---

## What the precomputed files contain

| File | Description |
|------|-------------|
| `stage1_cleaned_sorfs.csv` | Deduplicated sORF table from GRCh38 transcriptome scan |
| `*_sorted.bam` | STAR-aligned RNA-seq (Trimmomatic + STAR, GRCh38) |
| `*_psites.bed` | P-site BED with empirically calibrated offsets (chr-prefixed) |
