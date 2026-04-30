# Test Dataset — HeLa Asynchronous (Validation Run)

Runs the full pipeline on public HeLa asynchronous data and verifies
the expected output. Raw FASTQ files are downloaded automatically from
SRA — no manual file preparation needed.

**Runtime: ~2-4 hours** (dominated by SRA download + STAR alignment)

---

## Dataset

| Sample | SRR | Type |
|--------|-----|------|
| HeLa_async_RNA_rep1 | SRR3306581 | RNA-seq |
| HeLa_async_RNA_rep2 | SRR3306582 | RNA-seq |
| HeLa_async_RIBO_rep1 | SRR3306588 | Ribo-seq |
| HeLa_async_RIBO_rep2 | SRR3306589 | Ribo-seq |

Source: Aviner et al. (2017), *Nature Structural & Molecular Biology*, GSE79664.

---

## Prerequisites

**1. Install environment**
```bash
conda env create -f environment.yml
conda activate sorf-tool
```

**2. Reference files** — place in `data/raw/`:
```bash
# Genome FASTA
wget -P data/raw/ https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip data/raw/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# GTF
wget -P data/raw/ https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz
gunzip data/raw/Homo_sapiens.GRCh38.109.gtf.gz

# SortMeRNA rRNA databases
wget -P data/raw/rRNA_db/ https://github.com/sortmerna/sortmerna/raw/master/data/rRNA_databases/silva-euk-28s-id98.fasta
wget -P data/raw/rRNA_db/ https://github.com/sortmerna/sortmerna/raw/master/data/rRNA_databases/silva-euk-18s-id95.fasta
wget -P data/raw/rRNA_db/ https://github.com/sortmerna/sortmerna/raw/master/data/rRNA_databases/rfam-5s-database-id98.fasta
wget -P data/raw/rRNA_db/ https://github.com/sortmerna/sortmerna/raw/master/data/rRNA_databases/rfam-5.8s-database-id98.fasta
```

**3. Precomputed Phase 1 output** — download from [GitHub Releases](../../releases):
```
results/phase1/stage1_cleaned_sorfs.csv
```

---

## Run

```bash
bash test/run_test.sh
```

The script will:
1. Download FASTQs from SRA automatically (SRR3306581/82/88/89)
2. Run Trimmomatic trimming
3. Run STAR alignment (RNA-seq and Ribo-seq)
4. Run rRNA removal (Ribo-seq)
5. Calibrate P-site offsets empirically
6. Compute translation evidence and TE

---

## Expected output

`results/phase4/ribo_HeLa_async_common_translated_orfs.csv`

- **12 concordant sORFs** reproducible across both replicates
- **RPL26P19 (ORF_13030)** present — pseudogene-derived microprotein
  candidate reported in the paper

The test script checks both conditions and prints PASS/FAIL automatically.
