"""
Tier 3 (step 6) — validate on REAL AoU data at published ancestry-specific loci.

No ground truth here, so the test is reproducibility of KNOWN biology:
does the AMR-specific Tractor effect at an Amerindigenous-enriched locus move
toward its published direction/magnitude, with more power on the AMR track, when
the LAI came from the MXB panel?  And — the specificity check — do the AFR/EUR
control loci (DARC, SLC24A5, APOL1) NOT change with the AMR panel?

Runs against the real Tractor GWAS outputs produced for each panel (via the main
flow steps 9-12 + Tractor) on the real cohort + real phenotype.
"""

import argparse
import os
import pandas as pd
import numpy as np

# Parsed from eval_config.sh EVAL_LOCI (name|chrom|pos|trait|ancestry|role).
LOCI = [l.split("|") for l in os.environ.get("EVAL_LOCI", "").split("\n") if "|" in l]


def nearest_row(df, chrom, pos):
    sub = df[df["CHROM"].astype(str).str.replace("chr", "") == chrom.replace("chr", "")]
    if sub.empty:
        return None
    return sub.iloc[(sub["POS"] - pos).abs().argmin()]


def effect_at(gwas_path, chrom, pos, ancestry):
    df = pd.read_csv(gwas_path, sep="\t")
    df.columns = [c.strip() for c in df.columns]
    r = nearest_row(df, chrom, pos)
    if r is None:
        return {}
    return {
        f"{ancestry}_beta": r.get(f"{ancestry}_beta", np.nan),
        f"{ancestry}_neglog10p": -np.log10(max(r.get(f"{ancestry}_pval", 1), 1e-300)),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gwas", action="append", required=True,
                    help="label=path real-cohort Tractor GWAS (repeat per panel arm)")
    ap.add_argument("--loci", nargs="*", help="override: name,chrom,pos,ancestry per token")
    ap.add_argument("--out", default="realdata_positive_controls.tsv")
    args = ap.parse_args()

    panels = [spec.split("=", 1) for spec in args.gwas]     # [(label, path), ...]
    loci = LOCI if not args.loci else [t.split(",") for t in args.loci]
    rows = []
    for parts in loci:
        name, chrom, pos, ancestry = parts[0], parts[1], int(parts[2]), parts[4] if len(parts) > 4 else parts[3]
        row = {"locus": name, "chrom": chrom, "pos": pos, "expected_ancestry": ancestry}
        for label, path in panels:
            row.update({f"{label}_{k}": v for k, v in effect_at(path, chrom, pos, ancestry).items()})
        rows.append(row)

    out = pd.DataFrame(rows)
    out.to_csv(args.out, sep="\t", index=False)
    pd.set_option("display.width", 240)
    print(out.to_string(index=False))
    print("\nExpectation: at AMR loci (SLC16A11, ABCA1) the mxb_expanded panel gives "
          "a larger AMR-track -log10p and effect nearer the published value than the "
          "aou_default_1kg baseline; at AFR/EUR control loci all panels should agree "
          "(specificity — the AMR panel shouldn't move non-AMR signals).")


if __name__ == "__main__":
    main()
