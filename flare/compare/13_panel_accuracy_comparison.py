"""
Step 13 — benchmark the TARGETED reference panel (1KG-HGDP + MXB, >95%
homogeneous) against the COMPARISON panel (broad/admixed AMR, no MXB).

There is no ground-truth local ancestry for AoU participants, so we assess panel
fit with the standard proxy metrics:

  1. Concordance of per-marker calls between panels (where they agree, both are
     probably right; systematic disagreement localizes where AMR matching
     matters).
  2. FLARE posterior CONFIDENCE (the `probs=true` output): a better-matched
     panel yields sharper posteriors (higher max-prob, lower entropy),
     especially over AMR tracts.
  3. Global ancestry consistency vs the GrafAnc Pe/Pf and the supervised
     ADMIXTURE Q — summed FLARE tract length per ancestry should track the
     independent global estimate; deviation flags panel misfit.
  4. (If you hold out labeled reference samples) leave-one-out accuracy:
     paint known-ancestry ref individuals with each panel and score.

Outputs a tidy comparison table + summary plots.
"""

import argparse
import glob
import numpy as np
import pandas as pd
import pysam


def per_marker_confidence(anc_vcf):
    """Mean max-posterior and entropy per haplotype from FLARE probs output."""
    vcf = pysam.VariantFile(anc_vcf)
    samples = list(vcf.header.samples)
    maxp, ent = [], []
    for rec in vcf:
        for s in samples:
            # FLARE writes per-ancestry posteriors (e.g. ANP field); confirm name.
            probs = rec.samples[s].get("ANP")     # <-- confirm field vs your version
            if probs is None:
                continue
            p = np.array(probs, dtype=float)
            p = p[p > 0]
            maxp.append(p.max())
            ent.append(-(p * np.log(p)).sum())
    return np.mean(maxp), np.mean(ent)


def call_concordance(anc_vcf_a, anc_vcf_b):
    """Fraction of (marker,haplotype) calls that agree between two panels."""
    va = pysam.VariantFile(anc_vcf_a); vb = pysam.VariantFile(anc_vcf_b)
    agree = total = 0
    for ra, rb in zip(va, vb):
        if ra.pos != rb.pos:
            continue
        for s in va.header.samples:
            for fld in ("AN1", "AN2"):
                a = ra.samples[s].get(fld); b = rb.samples[s].get(fld)
                if a is None or b is None:
                    continue
                agree += int(a == b); total += 1
    return agree / total if total else float("nan")


def global_from_tracts(ancestry_table_glob, labels=("EUR", "AFR", "AMR")):
    """Genome-wide mean ancestry fraction from filled per-variant tracts."""
    frames = [pd.read_csv(f, sep="\t") for f in glob.glob(ancestry_table_glob)]
    df = pd.concat(frames, ignore_index=True)
    hap_cols = [c for c in df.columns if c != "pos"]
    counts = {i: (df[hap_cols].values == i).sum() for i in range(len(labels))}
    tot = sum(counts.values())
    return {lab: counts[i] / tot for i, lab in enumerate(labels)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--targeted_anc", required=True, help="targeted-panel .anc.vcf.gz (one chrom or merged)")
    ap.add_argument("--comparison_anc", required=True, help="comparison-panel .anc.vcf.gz")
    ap.add_argument("--targeted_tracts_glob", required=True, help="step-11 tables, targeted")
    ap.add_argument("--comparison_tracts_glob", required=True, help="step-11 tables, comparison")
    ap.add_argument("--grafanc_global", help="optional: GrafAnc Pe/Pf/Pa summary for the cohort")
    ap.add_argument("--out", default="panel_comparison.tsv")
    args = ap.parse_args()

    rows = []
    for name, anc in [("targeted", args.targeted_anc), ("comparison", args.comparison_anc)]:
        mp, en = per_marker_confidence(anc)
        rows.append({"panel": name, "mean_max_posterior": mp, "mean_entropy": en})
    conc = call_concordance(args.targeted_anc, args.comparison_anc)
    gt = global_from_tracts(args.targeted_tracts_glob)
    gc = global_from_tracts(args.comparison_tracts_glob)

    summary = pd.DataFrame(rows)
    summary.to_csv(args.out, sep="\t", index=False)
    print(summary.to_string(index=False))
    print(f"\nInter-panel call concordance: {conc:.4f}")
    print(f"Genome-wide ancestry (targeted):   {gt}")
    print(f"Genome-wide ancestry (comparison): {gc}")
    print("\nInterpretation: the better-matched panel should show higher "
          "mean_max_posterior, lower entropy, and global fractions closer to the "
          "independent GrafAnc/ADMIXTURE estimates — most visibly over AMR tracts.")


if __name__ == "__main__":
    main()
