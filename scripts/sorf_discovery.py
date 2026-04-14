#!/usr/bin/env python3
"""
sorf_discovery.py — Phase 1a: Scan transcriptome for sORFs.

Usage:
    python scripts/sorf_discovery.py \
        --fasta  data/raw/Homo_sapiens.GRCh38.cdna.all.fa \
        --gtf    data/raw/Homo_sapiens.GRCh38.115.gtf \
        --min_aa 10 \
        --max_aa 100 \
        --output results/phase1/stage1_novel_sorfs.csv
"""

import os
import re
import argparse
from typing import List, Tuple
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd

# ── CLI ──────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description="sORF discovery from transcriptome")
    p.add_argument("--fasta",   required=True, help="cDNA FASTA file")
    p.add_argument("--gtf",     required=True, help="Genome annotation GTF")
    p.add_argument("--min_aa",  type=int, default=10, help="Min peptide length (aa)")
    p.add_argument("--max_aa",  type=int, default=100, help="Max peptide length (aa)")
    p.add_argument("--output",  required=True, help="Output CSV path")
    return p.parse_args()

# ── DATA STRUCTURES ──────────────────────────────────────────────────────────

class Transcript:
    def __init__(self, chrom, strand):
        self.chrom  = chrom
        self.strand = strand
        self.exons  = []
        self.cds    = []

# ── PARSE GTF ────────────────────────────────────────────────────────────────

def parse_gtf(gtf_file: str):
    print(f"Parsing GTF: {gtf_file}")
    transcripts = {}
    with open(gtf_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue
            chrom, feature = fields[0], fields[2]
            start, end, strand = int(fields[3]), int(fields[4]), fields[6]
            m = re.search(r'transcript_id "([^"]+)"', fields[8])
            if not m:
                continue
            tid = m.group(1)
            if tid not in transcripts:
                transcripts[tid] = Transcript(chrom, strand)
            if feature == "exon":
                transcripts[tid].exons.append((start, end))
            if feature == "CDS":
                transcripts[tid].cds.append((start, end))
    for t in transcripts.values():
        t.exons.sort(key=lambda x: x[0])
    print(f"  Loaded {len(transcripts):,} transcripts")
    return transcripts

# ── LOAD FASTA ────────────────────────────────────────────────────────────────

def load_fasta(fasta_file: str):
    print(f"Loading FASTA: {fasta_file}")
    seqs = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        tid = record.id.split("|")[0]
        seqs[tid] = str(record.seq)
    print(f"  Loaded {len(seqs):,} sequences")
    return seqs

# ── COORDINATE MAPPING ────────────────────────────────────────────────────────

def transcript_to_genome(tid: str, t_start: int, t_end: int,
                         transcripts: dict) -> List[Tuple[str, int, int]]:
    t       = transcripts[tid]
    strand  = t.strand
    segments = []
    transcript_cursor = 0
    for exon_start, exon_end in t.exons:
        exon_length = exon_end - exon_start + 1
        exon_t_start = transcript_cursor
        exon_t_end   = transcript_cursor + exon_length
        if t_end <= exon_t_start:
            break
        if t_start >= exon_t_end:
            transcript_cursor += exon_length
            continue
        overlap_start = max(t_start, exon_t_start)
        overlap_end   = min(t_end,   exon_t_end)
        offset_start  = overlap_start - exon_t_start
        offset_end    = overlap_end   - exon_t_start
        if strand == "+":
            g_start = exon_start + offset_start
            g_end   = exon_start + offset_end - 1
        else:
            g_end   = exon_end - offset_start
            g_start = exon_end - offset_end + 1
        segments.append((t.chrom, g_start, g_end))
        transcript_cursor += exon_length
    return segments

# ── sORF SCANNING ─────────────────────────────────────────────────────────────

STOP_CODONS = {"TAA", "TAG", "TGA"}

def scan_sorfs(fasta_sequences: dict, transcripts: dict,
               min_aa: int, max_aa: int):
    print(f"Scanning sORFs (min_aa={min_aa}, max_aa={max_aa})...")
    results = []
    total   = len(fasta_sequences)
    for idx, (tid, seq) in enumerate(fasta_sequences.items(), 1):
        if idx % 5000 == 0:
            print(f"  {idx:,}/{total:,} transcripts | sORFs: {len(results):,}")
        if tid not in transcripts:
            continue
        if len(seq) < min_aa * 3:
            continue
        if "ATG" not in seq:
            continue
        seq      = seq.upper()
        seq_len  = len(seq)
        t_data   = transcripts[tid]
        cds_regs = t_data.cds
        for frame in range(3):
            i = frame
            while i < seq_len - 2:
                if seq[i:i+3] != "ATG":
                    i += 3
                    continue
                j = i + 3
                while j < seq_len - 2:
                    stop = seq[j:j+3]
                    if stop in STOP_CODONS:
                        aa_len = (j + 3 - i) // 3 - 1
                        if min_aa <= aa_len <= max_aa:
                            nt_seq   = seq[i:j+3]
                            segments = transcript_to_genome(tid, i, j+3, transcripts)
                            overlaps = any(
                                g_start <= cds_end and g_end >= cds_start
                                for cds_start, cds_end in cds_regs
                                for _, g_start, g_end in segments
                            )
                            if not overlaps:
                                results.append({
                                    "transcript_id":    tid,
                                    "chromosome":       t_data.chrom,
                                    "strand":           t_data.strand,
                                    "transcript_start": i,
                                    "transcript_end":   j + 3,
                                    "aa_length":        aa_len,
                                    "genomic_segments": segments,
                                    "nt_sequence":      nt_seq,
                                    "aa_sequence":      str(Seq(nt_seq).translate(to_stop=True)),
                                })
                        i = j + 3
                        break
                    j += 3
                else:
                    i += 3
    print(f"  Total sORFs found: {len(results):,}")
    return results

# ── MAIN ──────────────────────────────────────────────────────────────────────

def main():
    args        = parse_args()
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    transcripts = parse_gtf(args.gtf)
    seqs        = load_fasta(args.fasta)
    results     = scan_sorfs(seqs, transcripts, args.min_aa, args.max_aa)
    df          = pd.DataFrame(results)
    df.to_csv(args.output, index=False)
    print(f"\nSaved {len(df):,} sORFs → {args.output}")

if __name__ == "__main__":
    main()
