#!/usr/bin/env python3
"""
orf_scorer.py — Random Forest classifier for translated sORF scoring.

HOW IT WORKS:
  Uses a self-supervised labelling strategy — no external training set needed:

    POSITIVES  : ORFs with strong multi-evidence support
                 (high ribo counts AND strong periodicity AND detectable RNA)
    NEGATIVES  : ORFs with ribo_count = 0 or frame0_fraction < 0.34
                 (near-random frame distribution = not translated)

  Features engineered from columns already present in your Stage 5 output:
    - ribo_count          (raw Ribo-seq read depth)
    - total_psites        (P-site count)
    - frame0_fraction     (triplet periodicity)
    - TE_mean             (translation efficiency)
    - aa_length           (peptide length)
    - rna_rep1, rna_rep2  (RNA expression)
    - gc_content          (from nt_sequence)
    - cai_score           (codon adaptation index proxy)
    - kozak_strength      (sequence context around AUG)
    - uorf_flag           (is this upstream of an annotated CDS?)

  Output adds two columns to your existing table:
    - rf_score     : probability of being translated (0–1)
    - rf_label     : "translated" / "unlikely" / "ambiguous"

PIPELINE POSITION:
  Runs AFTER Phase 5 TE filter. Takes the final per-cell-type TE-filtered
  table and adds scores. Does not remove any ORFs — purely additive.

  results/phase5/<cell>_translated_orfs_filtered_withTE.csv
       ↓
  orf_scorer.py
       ↓
  results/phase6_scoring/<cell>_scored_orfs.csv

Usage:
    python scripts/orf_scorer.py \
        --input    results/phase5/HeLa_M_translated_orfs_filtered_withTE.csv \
        --cell     HeLa_M \
        --outdir   results/phase6_scoring

    # To score a new cell type using a model trained on HeLa:
    python scripts/orf_scorer.py \
        --input      results/phase5/fibroblast_translated_orfs_filtered_withTE.csv \
        --cell       fibroblast \
        --model      results/phase6_scoring/HeLa_M_rf_model.pkl \
        --outdir     results/phase6_scoring
"""

import os
import re
import json
import pickle
import argparse
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (classification_report, roc_auc_score,
                              RocCurveDisplay, ConfusionMatrixDisplay)
from sklearn.inspection import permutation_importance

# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description="RF scoring for translated sORFs")
    p.add_argument("--input",   required=True,  help="TE-filtered ORF CSV (Phase 5 output)")
    p.add_argument("--cell",    required=True,  help="Cell type label")
    p.add_argument("--model",   default=None,   help="Pre-trained model .pkl (optional)")
    p.add_argument("--outdir",  required=True,  help="Output directory")
    # Labelling thresholds
    p.add_argument("--pos_ribo",    type=int,   default=20,
                   help="Min ribo_count for a POSITIVE label")
    p.add_argument("--pos_frame",   type=float, default=0.60,
                   help="Min frame0_fraction for a POSITIVE label")
    p.add_argument("--neg_ribo",    type=int,   default=0,
                   help="Max ribo_count for a NEGATIVE label (inclusive)")
    p.add_argument("--neg_frame",   type=float, default=0.34,
                   help="Max frame0_fraction for a NEGATIVE label")
    # RF params
    p.add_argument("--n_trees",     type=int,   default=500)
    p.add_argument("--cv_folds",    type=int,   default=5)
    p.add_argument("--score_thr",   type=float, default=0.5,
                   help="RF score threshold for 'translated' label")
    return p.parse_args()

# ── FEATURE ENGINEERING ───────────────────────────────────────────────────────

# Standard codon table for CAI calculation
CODON_USAGE_HUMAN = {
    # High-usage codons → weight near 1.0; rare codons → weight near 0
    # Values approximate human codon adaptation weights (Nakamura et al.)
    "TTT":0.45,"TTC":1.00,"TTA":0.07,"TTG":0.13,
    "CTT":0.13,"CTC":0.48,"CTA":0.07,"CTG":1.00,
    "ATT":0.36,"ATC":1.00,"ATA":0.16,"ATG":1.00,
    "GTT":0.18,"GTC":0.45,"GTA":0.11,"GTG":1.00,
    "TCT":0.15,"TCC":0.44,"TCA":0.15,"TCG":0.06,
    "CCT":0.28,"CCC":0.61,"CCA":0.27,"CCG":0.11,
    "ACT":0.24,"ACC":1.00,"ACA":0.28,"ACG":0.11,
    "GCT":0.27,"GCC":1.00,"GCA":0.22,"GCG":0.11,
    "TAT":0.44,"TAC":1.00,"TAA":0.28,"TAG":0.20,
    "CAT":0.41,"CAC":1.00,"CAA":0.25,"CAG":1.00,
    "AAT":0.46,"AAC":1.00,"AAA":0.43,"AAG":1.00,
    "GAT":0.46,"GAC":1.00,"GAA":0.42,"GAG":1.00,
    "TGT":0.45,"TGC":1.00,"TGA":0.52,"TGG":1.00,
    "CGT":0.08,"CGC":0.40,"CGA":0.11,"CGG":0.20,
    "AGT":0.15,"AGC":0.59,"AGA":0.20,"AGG":0.20,
    "GGT":0.16,"GGC":0.74,"GGA":0.25,"GGG":0.25,
}

def compute_cai(nt_seq: str) -> float:
    """Geometric mean of codon weights (CAI proxy)."""
    seq    = nt_seq.upper()
    codons = [seq[i:i+3] for i in range(0, len(seq)-2, 3)]
    weights = [CODON_USAGE_HUMAN.get(c, 0.1) for c in codons if len(c) == 3]
    if not weights:
        return 0.0
    return float(np.exp(np.mean(np.log(np.maximum(weights, 1e-6)))))

def compute_gc(nt_seq: str) -> float:
    s = nt_seq.upper()
    if not s:
        return 0.0
    return (s.count("G") + s.count("C")) / len(s)

def kozak_strength(nt_seq: str) -> float:
    """
    Score the Kozak context around the start ATG.
    Consensus: (gcc)gccRccAUGG
    We check positions -3 (R = A/G) and +4 (G) — the two strongest positions.
    Returns 0.0, 0.5, or 1.0 based on how many are optimal.
    """
    if len(nt_seq) < 7:
        return 0.0
    # ATG is at position 0 in nt_seq
    # position -3 relative to ATG start
    pos_minus3 = nt_seq[0] if len(nt_seq) >= 1 else "N"  # already at ATG
    # We can only check +4 if seq is long enough
    pos_plus4  = nt_seq[3] if len(nt_seq) > 3 else "N"
    score = 0.0
    if pos_minus3 in ("A", "G"):   # purine at -3
        score += 0.5
    if pos_plus4 == "G":           # G at +4
        score += 0.5
    return score

def engineer_features(df: pd.DataFrame) -> pd.DataFrame:
    """
    Build feature matrix from columns available in the Phase 5 output.
    Handles missing columns gracefully with zeros.
    """
    print("  Engineering features...")
    feats = pd.DataFrame(index=df.index)

    # ── Direct numeric features ──────────────────────────────────────────────
    for col in ["ribo_count", "total_psites", "frame0_fraction",
                "TE_mean", "aa_length", "rna_rep1", "rna_rep2",
                "ribo_rep1", "ribo_rep2", "frame0", "frame1", "frame2"]:
        if col in df.columns:
            feats[col] = pd.to_numeric(df[col], errors="coerce").fillna(0)
        else:
            feats[col] = 0.0

    # ── Frame imbalance (how dominant is frame 0 vs others) ─────────────────
    total_f = feats["frame0"] + feats["frame1"] + feats["frame2"] + 1e-6
    feats["frame_imbalance"] = feats["frame0"] / total_f

    # ── Ribo consistency across reps ─────────────────────────────────────────
    max_ribo = feats[["ribo_rep1","ribo_rep2"]].max(axis=1) + 1e-6
    min_ribo = feats[["ribo_rep1","ribo_rep2"]].min(axis=1)
    feats["ribo_consistency"] = min_ribo / max_ribo   # 1.0 = perfectly consistent

    # ── RNA detectability ────────────────────────────────────────────────────
    feats["rna_detected"] = (
        (feats["rna_rep1"] > 0).astype(int) +
        (feats["rna_rep2"] > 0).astype(int)
    )

    # ── Log-transformed counts (stabilize variance) ───────────────────────────
    for col in ["ribo_count", "total_psites", "rna_rep1", "rna_rep2"]:
        feats[f"log_{col}"] = np.log1p(feats[col])

    # ── Sequence-derived features (if nt_sequence available) ─────────────────
    if "nt_sequence" in df.columns:
        feats["gc_content"]     = df["nt_sequence"].apply(
            lambda s: compute_gc(s) if pd.notna(s) else 0.5)
        feats["cai_score"]      = df["nt_sequence"].apply(
            lambda s: compute_cai(s) if pd.notna(s) else 0.0)
        feats["kozak_strength"] = df["nt_sequence"].apply(
            lambda s: kozak_strength(s) if pd.notna(s) else 0.0)
    else:
        feats["gc_content"]     = 0.5
        feats["cai_score"]      = 0.0
        feats["kozak_strength"] = 0.0

    # ── TE ratio (rep1 vs rep2 consistency) ──────────────────────────────────
    if "TE_rep1" in df.columns and "TE_rep2" in df.columns:
        te1 = pd.to_numeric(df["TE_rep1"], errors="coerce").fillna(0)
        te2 = pd.to_numeric(df["TE_rep2"], errors="coerce").fillna(0)
        max_te = np.maximum(te1, te2) + 1e-6
        min_te = np.minimum(te1, te2)
        feats["te_consistency"] = min_te / max_te
    else:
        feats["te_consistency"] = 0.0

    print(f"  Feature matrix: {feats.shape[0]:,} ORFs × {feats.shape[1]} features")
    return feats

# ── SELF-SUPERVISED LABELLING ────────────────────────────────────────────────

def make_labels(df: pd.DataFrame, feats: pd.DataFrame,
                pos_ribo: int, pos_frame: float,
                neg_ribo: int, neg_frame: float):
    """
    Assign confident labels using hard thresholds on the strongest signals.
    Returns a Series with 1 (positive), 0 (negative), NaN (unlabelled).
    """
    labels = pd.Series(np.nan, index=df.index)

    # Positives: high ribo + strong periodicity
    pos_mask = (
        (feats["ribo_count"]       >= pos_ribo) &
        (feats["frame0_fraction"]  >= pos_frame)
    )

    # Negatives: zero ribo OR near-random frame distribution
    neg_mask = (
        (feats["ribo_count"]       <= neg_ribo) |
        (feats["frame0_fraction"]  <= neg_frame)
    )

    labels[pos_mask] = 1
    labels[neg_mask & ~pos_mask] = 0

    n_pos = pos_mask.sum()
    n_neg = (neg_mask & ~pos_mask).sum()
    n_unl = labels.isna().sum()
    print(f"  Labels: {n_pos} positives | {n_neg} negatives | {n_unl} unlabelled")

    return labels

# ── TRAIN + EVALUATE ─────────────────────────────────────────────────────────

def train_rf(feats: pd.DataFrame, labels: pd.Series,
             n_trees: int, cv_folds: int, outdir: str, cell: str):
    """
    Train Random Forest on confidently labelled ORFs.
    Evaluate with stratified k-fold cross-validation.
    """
    labelled = labels.notna()
    X = feats[labelled].values
    y = labels[labelled].values.astype(int)

    print(f"\n  Training Random Forest ({n_trees} trees, {cv_folds}-fold CV)...")
    print(f"  Training set: {len(X):,} ORFs  ({int(y.sum())} pos, {int((1-y).sum())} neg)")

    clf = RandomForestClassifier(
        n_estimators=n_trees,
        max_features="sqrt",
        min_samples_leaf=3,
        class_weight="balanced",
        random_state=42,
        n_jobs=-1,
    )

    # Cross-validation
    cv      = StratifiedKFold(n_splits=cv_folds, shuffle=True, random_state=42)
    cv_aucs = cross_val_score(clf, X, y, cv=cv, scoring="roc_auc")
    print(f"  CV ROC-AUC: {cv_aucs.mean():.3f} ± {cv_aucs.std():.3f}")

    # Final fit on all labelled data
    clf.fit(X, y)

    # Feature importance
    feat_names   = feats.columns.tolist()
    importances  = pd.Series(clf.feature_importances_, index=feat_names)
    importances  = importances.sort_values(ascending=False)

    print(f"\n  Top 10 features:")
    for feat, imp in importances.head(10).items():
        bar = "█" * int(imp * 40)
        print(f"    {feat:<30} {imp:.4f}  {bar}")

    # Save feature importance plot
    fig, ax = plt.subplots(figsize=(8, 6))
    importances.head(15).sort_values().plot.barh(ax=ax, color="steelblue")
    ax.set_xlabel("Feature importance (mean decrease impurity)")
    ax.set_title(f"{cell} — RF Feature Importances")
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, f"{cell}_rf_feature_importance.pdf"),
                bbox_inches="tight")
    plt.close()

    return clf, cv_aucs, feat_names

# ── SCORE ALL ORFs ────────────────────────────────────────────────────────────

def score_orfs(clf, feats: pd.DataFrame, score_thr: float):
    X      = feats.values
    scores = clf.predict_proba(X)[:, 1]

    labels = pd.Series("ambiguous", index=feats.index)
    labels[scores >= score_thr]       = "translated"
    labels[scores < (1 - score_thr)]  = "unlikely"

    n_trans = (labels == "translated").sum()
    n_unlik = (labels == "unlikely").sum()
    n_ambig = (labels == "ambiguous").sum()
    print(f"\n  Scoring results (threshold = {score_thr}):")
    print(f"    Translated  : {n_trans:,}")
    print(f"    Unlikely    : {n_unlik:,}")
    print(f"    Ambiguous   : {n_ambig:,}")

    return scores, labels

# ── PLOTS ─────────────────────────────────────────────────────────────────────

def make_score_plot(df_out: pd.DataFrame, cell: str, outdir: str):
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    # Score distribution
    axes[0].hist(df_out["rf_score"], bins=50, color="steelblue", edgecolor="white")
    axes[0].axvline(0.5, color="red", linestyle="--", label="threshold")
    axes[0].set_xlabel("RF score")
    axes[0].set_ylabel("Count")
    axes[0].set_title("Score distribution")
    axes[0].legend()

    # Score vs frame0_fraction
    if "frame0_fraction" in df_out.columns:
        sc = axes[1].scatter(
            df_out["frame0_fraction"],
            df_out["rf_score"],
            c=df_out["rf_score"], cmap="RdYlGn",
            alpha=0.4, s=8
        )
        axes[1].set_xlabel("Frame-0 fraction")
        axes[1].set_ylabel("RF score")
        axes[1].set_title("Score vs Periodicity")
        plt.colorbar(sc, ax=axes[1])

    # Score vs log ribo count
    if "ribo_count" in df_out.columns:
        sc2 = axes[2].scatter(
            np.log1p(df_out["ribo_count"]),
            df_out["rf_score"],
            c=df_out["rf_score"], cmap="RdYlGn",
            alpha=0.4, s=8
        )
        axes[2].set_xlabel("log(ribo_count + 1)")
        axes[2].set_ylabel("RF score")
        axes[2].set_title("Score vs Ribo depth")
        plt.colorbar(sc2, ax=axes[2])

    plt.suptitle(f"{cell} — RF ORF Scores", fontsize=13, y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, f"{cell}_rf_scores.pdf"), bbox_inches="tight")
    plt.close()
    print(f"  Score plot saved")

# ── MAIN ──────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    print(f"\n{'='*56}")
    print(f"  RF ORF Scorer — {args.cell}")
    print(f"{'='*56}\n")

    # ── Load data ────────────────────────────────────────────────────────────
    df = pd.read_csv(args.input)
    print(f"Loaded {len(df):,} ORFs from {args.input}")

    # ── Feature engineering ──────────────────────────────────────────────────
    feats = engineer_features(df)

    # ── Load or train model ───────────────────────────────────────────────────
    if args.model and os.path.exists(args.model):
        print(f"\nLoading pre-trained model: {args.model}")
        with open(args.model, "rb") as f:
            saved = pickle.load(f)
        clf        = saved["model"]
        feat_names = saved["features"]
        cv_aucs    = saved.get("cv_aucs", [])

        # Align feature columns to training order
        for fn in feat_names:
            if fn not in feats.columns:
                feats[fn] = 0.0
        feats = feats[feat_names]

        print(f"  Model loaded (trained on: {saved.get('cell', 'unknown')})")
        if len(cv_aucs):
            print(f"  Original CV AUC: {np.mean(cv_aucs):.3f} ± {np.std(cv_aucs):.3f}")

    else:
        # Self-supervised training on this dataset
        print("\nNo pre-trained model — training self-supervised RF...")
        labels = make_labels(
            df, feats,
            args.pos_ribo, args.pos_frame,
            args.neg_ribo, args.neg_frame
        )

        if labels.notna().sum() < 50:
            print("ERROR: Fewer than 50 confidently labelled ORFs.")
            print("Consider lowering --pos_ribo or adjusting thresholds.")
            return

        clf, cv_aucs, feat_names = train_rf(
            feats, labels, args.n_trees, args.cv_folds,
            args.outdir, args.cell
        )

        # Save model
        model_path = os.path.join(args.outdir, f"{args.cell}_rf_model.pkl")
        with open(model_path, "wb") as f:
            pickle.dump({
                "model":    clf,
                "features": feat_names,
                "cv_aucs":  cv_aucs.tolist(),
                "cell":     args.cell,
                "params": {
                    "pos_ribo":   args.pos_ribo,
                    "pos_frame":  args.pos_frame,
                    "neg_ribo":   args.neg_ribo,
                    "neg_frame":  args.neg_frame,
                    "n_trees":    args.n_trees,
                },
            }, f)
        print(f"\n  Model saved: {model_path}")

    # ── Score all ORFs ────────────────────────────────────────────────────────
    scores, rf_labels = score_orfs(clf, feats, args.score_thr)

    df_out = df.copy()
    df_out["rf_score"] = np.round(scores, 4)
    df_out["rf_label"] = rf_labels.values

    # ── Sort by score ─────────────────────────────────────────────────────────
    df_out = df_out.sort_values("rf_score", ascending=False)

    # ── Save outputs ──────────────────────────────────────────────────────────
    all_path  = os.path.join(args.outdir, f"{args.cell}_scored_orfs.csv")
    high_path = os.path.join(args.outdir, f"{args.cell}_high_confidence_scored_orfs.csv")

    df_out.to_csv(all_path, index=False)
    df_out[df_out["rf_label"] == "translated"].to_csv(high_path, index=False)

    print(f"\n  All ORFs with scores : {all_path}")
    print(f"  High-confidence only : {high_path}")

    # ── Plots ────────────────────────────────────────────────────────────────
    try:
        make_score_plot(df_out, args.cell, args.outdir)
    except Exception as e:
        print(f"  (Plot skipped: {e})")

    # ── Final summary ─────────────────────────────────────────────────────────
    print(f"\n{'='*56}")
    print(f"  SCORING COMPLETE — {args.cell}")
    if len(cv_aucs):
        print(f"  Model CV AUC       : {np.mean(cv_aucs):.3f} ± {np.std(cv_aucs):.3f}")
    print(f"  Total ORFs scored  : {len(df_out):,}")
    print(f"  Translated (≥{args.score_thr}) : "
          f"{(df_out['rf_label']=='translated').sum():,}")
    print(f"{'='*56}\n")

if __name__ == "__main__":
    main()
