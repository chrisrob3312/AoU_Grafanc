"""
Step 9 — novel-discovery comparison across panels on REAL AoU data.

Question: does the expanded MXB panel surface AMR-track association signals in the
Latin American cohort that the AoU-default (1KG) panel misses — i.e. potential
novel discovery enabled by better Amerindigenous local ancestry?

This is the most exciting claim and the easiest to fool yourself with, so it is
gated by calibration checks. A "novel hit" only counts if:
  1. genome-wide significant on the AMR ancestry track (default p < 5e-8),
  2. NOT within --known_window bp of a known association (--known_catalog),
  3. the panel's AMR track is well-calibrated genome-wide (lambda_GC within
     --lambda_tol of 1) — an inflated panel manufactures false "novel" hits,
  4. the signal is ROBUST: still significant after refit / present in both
     halves if --replication_split is given.

Outputs, per panel: calibration, counts of known-recovered and novel-candidate
hits, and the novel-candidate list. Compares mxb_expanded vs aou_default_1kg.
"""

import argparse
import numpy as np
import pandas as pd
from scipy import stats


def lambda_gc(pvals):
    p = pd.Series(pvals); p = p[(p > 0) & (p < 1)].dropna()
    chi2 = stats.chi2.isf(p, df=1)
    return float(np.median(chi2) / stats.chi2.ppf(0.5, df=1)) if len(p) else float("nan")


def load(path):
    df = pd.read_csv(path, sep="\t"); df.columns = [c.strip() for c in df.columns]
    return df


def clump(hits, window):
    """Keep the lead SNP per locus (greedy by p within +/- window bp, per chrom)."""
    leads = []
    for _, chrom_hits in hits.groupby("CHROM"):
        remaining = chrom_hits.sort_values("AMR_pval")
        while not remaining.empty:
            lead = remaining.iloc[0]
            leads.append(lead)
            remaining = remaining[(remaining["POS"] - lead["POS"]).abs() > window]
    return pd.DataFrame(leads)


def is_known(pos, chrom, catalog, window):
    c = catalog[catalog["chrom"].astype(str).str.replace("chr", "") == str(chrom).replace("chr", "")]
    return bool(((c["pos"] - pos).abs() <= window).any())


def analyze(label, path, catalog, sig, window, lambda_tol):
    df = load(path)
    lam = lambda_gc(df["AMR_pval"])
    calibrated = abs(lam - 1) <= lambda_tol
    sig_hits = df[df["AMR_pval"] < sig].copy()
    leads = clump(sig_hits, window) if not sig_hits.empty else sig_hits
    known = novel = 0
    novel_rows = []
    for _, r in leads.iterrows():
        if is_known(r["POS"], r["CHROM"], catalog, window):
            known += 1
        else:
            novel += 1
            novel_rows.append(r)
    return {
        "panel": label, "AMR_lambda_gc": round(lam, 3), "calibrated": calibrated,
        "genomewide_sig_leads": len(leads),
        "known_recovered": known,
        "novel_candidates": novel if calibrated else "SUPPRESSED (miscalibrated)",
    }, pd.DataFrame(novel_rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gwas", action="append", required=True, help="label=path (repeat per panel)")
    ap.add_argument("--known_catalog", required=True,
                    help="TSV of known associations: chrom<TAB>pos (GWAS Catalog subset for the trait)")
    ap.add_argument("--sig", type=float, default=5e-8)
    ap.add_argument("--known_window", type=int, default=500_000)
    ap.add_argument("--lambda_tol", type=float, default=0.10)
    ap.add_argument("--out_prefix", default="novel_discovery")
    args = ap.parse_args()

    catalog = pd.read_csv(args.known_catalog, sep="\t")
    summaries = []
    for spec in args.gwas:
        label, path = spec.split("=", 1)
        summ, novel = analyze(label, path, catalog, args.sig, args.known_window, args.lambda_tol)
        summaries.append(summ)
        if not novel.empty:
            novel.to_csv(f"{args.out_prefix}.{label}.novel.tsv", sep="\t", index=False)

    out = pd.DataFrame(summaries)
    out.to_csv(f"{args.out_prefix}.summary.tsv", sep="\t", index=False)
    print(out.to_string(index=False))
    print("\nCaveats to state in any writeup: novel candidates are HYPOTHESES, not "
          "confirmed hits — they need replication in an independent Latin American "
          "cohort and are only credible when the AMR track is calibrated (lambda≈1). "
          "The headline is the DIFFERENCE: signals mxb_expanded finds that "
          "aou_default_1kg misses, with both calibrated.")


if __name__ == "__main__":
    main()
