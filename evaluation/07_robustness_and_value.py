"""
Tier 4 (step 7) — genome-wide robustness + broader value arguments for the
MXB-augmented AMR panel. These don't need injected truth; they characterize
behavior across the whole genome / cohort.

Reports:
  1. Genome-wide AMR-track lambda_GC per panel (calibration; mis-assigned tracts
     inflate/deflate the AMR track).
  2. Admixture-mapping peak: association of AMR LOCAL ANCESTRY dosage with the
     trait, at the primary AMR locus — sharper/better-localized with a better panel.
  3. Global-vs-local ancestry consistency: summed AMR local ancestry per person
     vs the independent GrafAnc/ADMIXTURE global AMR estimate (better panel →
     tighter agreement, less systematic AMR mis-calling).
  4. Tract-length distribution: excess short AMR tracts = fragmentation from panel
     misfit.

Each is a separate function so you can run only what you have inputs for.
"""

import argparse
import numpy as np
import pandas as pd
from scipy import stats


def lambda_gc(pvals):
    p = pd.Series(pvals)
    p = p[(p > 0) & (p < 1)].dropna()
    chi2 = stats.chi2.isf(p, df=1)
    return float(np.median(chi2) / stats.chi2.ppf(0.5, df=1)) if len(p) else float("nan")


def admixture_map(local_anc_amr_csv, pheno_tsv):
    """Regress trait on per-person AMR local-ancestry dosage at each position."""
    la = pd.read_csv(local_anc_amr_csv)          # pos + one col per person (AMR dosage 0-2)
    phe = pd.read_csv(pheno_tsv, sep="\t").set_index("person_id")["y"]
    people = [c for c in la.columns if c in phe.index]
    y = phe.loc[people].values
    out = []
    for _, r in la.iterrows():
        x = r[people].values.astype(float)
        if np.nanstd(x) == 0:
            continue
        b, _, rval, pval, _ = stats.linregress(x, y)
        out.append((r["pos"], b, pval))
    df = pd.DataFrame(out, columns=["pos", "beta", "pval"])
    return df


def global_local_consistency(local_anc_amr_csv, global_amr_tsv):
    """Mean AMR local ancestry per person vs their global AMR estimate."""
    la = pd.read_csv(local_anc_amr_csv)
    people = [c for c in la.columns if c != "pos"]
    local_mean = la[people].mean(axis=0) / 2.0     # dosage 0-2 -> fraction
    g = pd.read_csv(global_amr_tsv, sep="\t").set_index("person_id")["global_AMR"]
    common = [p for p in people if p in g.index]
    r = np.corrcoef(local_mean.loc[common], g.loc[common])[0, 1]
    rmse = float(np.sqrt(np.mean((local_mean.loc[common] - g.loc[common]) ** 2)))
    return {"local_vs_global_r": round(r, 4), "local_vs_global_rmse": round(rmse, 4)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", required=True, choices=["targeted", "comparison"])
    ap.add_argument("--gwas", help="Tractor GWAS for genome-wide AMR lambda")
    ap.add_argument("--local_anc_amr", help="per-variant AMR local-ancestry dosage csv")
    ap.add_argument("--pheno", help="phenotype tsv (for admixture mapping)")
    ap.add_argument("--global_amr", help="per-person global AMR estimate tsv")
    ap.add_argument("--out_prefix", default="robustness")
    args = ap.parse_args()

    summary = {"panel": args.panel}
    if args.gwas:
        g = pd.read_csv(args.gwas, sep="\t"); g.columns = [c.strip() for c in g.columns]
        summary["AMR_lambda_gc"] = round(lambda_gc(g["AMR_pval"]), 3)
    if args.local_anc_amr and args.pheno:
        am = admixture_map(args.local_anc_amr, args.pheno)
        am.to_csv(f"{args.out_prefix}.{args.panel}.admixmap.tsv", sep="\t", index=False)
        summary["admixmap_peak_neglog10p"] = round(-np.log10(max(am["pval"].min(), 1e-300)), 3)
    if args.local_anc_amr and args.global_amr:
        summary.update(global_local_consistency(args.local_anc_amr, args.global_amr))

    pd.DataFrame([summary]).to_csv(f"{args.out_prefix}.{args.panel}.summary.tsv", sep="\t", index=False)
    print(summary)


if __name__ == "__main__":
    main()
