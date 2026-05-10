# Test Dataset — HeLa Asynchronous (Validation Run)

Validates the pipeline using precomputed HeLa asynchronous data from Zenodo.
Only Phases 4 and 5 are run. No alignment or download of raw FASTQs needed.

**Zenodo DOI: https://doi.org/10.5281/zenodo.20084600**

**Runtime: ~10 minutes**

---

## Dataset

HeLa asynchronous cells, 2 biological replicates (RNA-seq + Ribo-seq).
Source: Aviner et al. (2017), *Nature Structural & Molecular Biology*, GSE79664.

Precomputed files provided:
- `stage1_cleaned_sorfs.csv` — sORF table from GRCh38 transcriptome scan
- `HeLa_async_RNA_rep1/2_sorted.bam` — STAR-aligned RNA-seq BAMs
- `HeLa_async_RIBO_rep1/2_psites.bed` — P-site BEDs with empirical offsets

---

## Setup

**Step 1 — Install environment:**
```bash
conda env create -f environment.yml
conda activate ribowin
```

**Step 2 — Download GTF:**
```bash
wget -P data/raw/ https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz
gunzip data/raw/Homo_sapiens.GRCh38.115.gtf.gz
```

**Step 3 — Run (Zenodo data downloads automatically):**
```bash
bash test/run_test.sh
```

The script downloads the Zenodo archive automatically if files are missing.

---

## Expected output

`results/phase4/ribo_HeLa_async_common_translated_orfs.csv`

- **12 concordant sORFs** reproducible across both replicates
- **RPL26P19** (chr5:56,504,xxx) present — pseudogene-derived microprotein

The test script validates both conditions and prints PASS/FAIL.
