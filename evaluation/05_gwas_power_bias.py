"""
Tier 2 (step 5) — compare the two panels' Tractor GWAS against the KNOWN truth.

Because step 3 injected the per-ancestry betas, we can measure not just power but
BIAS — the thing that really shows panel value. Metrics per panel:

  * power           : -log10 p on the AMR ancestry track at the causal variant
  * AMR beta bias   : (estimated AMR effect) - (injected AMR effect)
  * AMR beta RMSE   : across replicate seeds if provided
  * leakage         : |estimated EUR/AFR effect| at the causal site (should be ~0;
                      mis-assigned AMR tracts leak the effect into other tracks)
  * calibration     : genomic-control lambda of the AMR track genome-wide
  * localization    : distance (bp) from the AMR-track peak to the true causal pos

The MXB panel should show higher power, smaller AMR bias, smaller leakage, and
lambda closer to 1.
"""

import argparse
import glob
import numpy as np
import pandas as pd
from scipy import stats


def load_tractor(path):
    """Tractor output → DataFrame with per-ancestry effect + p columns.
    Column names vary by Tractor version; map them here."""
    df = pd.read_csv(path, sep="\t")
    df.columns = [c.strip() for c in df.columns]
    return df  # expects e.g. POS, AMR_beta, AMR_pval, EUR_beta, AFR_beta, ...


def lambda_gc(pvals):
    p = pvals[(pvals > 0) & (pvals < 1)].dropna()
    if p.empty:
        return float("nan")
    chi2 = stats.chi2.isf(p, df=1)
    return float(np.median(chi2) / stats.chi2.ppf(0.5, df=1))


def summarize(path, causal_pos, truth):
    df = load_tractor(path)
    row = df.iloc[(df["POS"] - causal_pos).abs().argmin()]
    amr_peak_pos = df.loc[df["AMR_pval"].idxmin(), "POS"]
    return {
        "AMR_power_neglog10p": round(-np.log10(max(row["AMR_pval"], 1e-300)), 3),
        "AMR_beta_est": round(row["AMR_beta"], 4),
        "AMR_beta_bias": round(row["AMR_beta"] - truth["beta_AMR"], 4),
        "EUR_leakage_abs": round(abs(row.get("EUR_beta", np.nan)), 4),
        "AFR_leakage_abs": round(abs(row.get("AFR_beta", np.nan)), 4),
        "AMR_lambda_gc": round(lambda_gc(df["AMR_pval"]), 3),
        "AMR_peak_dist_bp": int(abs(amr_peak_pos - causal_pos)),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gwas_targeted", required=True)
    ap.add_argument("--gwas_comparison", required=True)
    ap.add_argument("--truth", required=True, help="*.truth.tsv from step 3")
    ap.add_argument("--out", default="panel_gwas_comparison.tsv")
    args = ap.parse_args()

    truth = pd.read_csv(args.truth, sep="\t").iloc[0]
    causal = int(truth["causal_pos"])
    rows = []
    for name, path in [("targeted", args.gwas_targeted), ("comparison", args.gwas_comparison)]:
        rows.append({"panel": name, **summarize(path, causal, truth)})
    out = pd.DataFrame(rows)
    out.to_csv(args.out, sep="\t", index=False)
    print(f"Injected AMR beta = {truth['beta_AMR']}")
    print(out.to_string(index=False))
    print("\nBetter panel: higher AMR_power, smaller |AMR_beta_bias| and leakage, "
          "lambda≈1, smaller peak distance.")


if __name__ == "__main__":
    main()
