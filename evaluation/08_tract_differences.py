"""
Core comparison (step 8) — do we get DIFFERENT tracts across reference panels?

This directly answers "does the expanded MXB panel change the local-ancestry
tracts for a Latin American cohort, vs a small-homogeneous AMR panel, a
large-admixed AMR panel, and an AoU-default 1KG panel?" It works on REAL AoU
FLARE output (no ground truth needed) — the point is to show the panel choice
materially changes the painting, and where.

For every pair of panels it reports, over aligned (position, haplotype) calls:
  * overall discordance fraction,
  * AMR-focused discordance: how often one panel calls AMR where the other does
    not (the ancestry that the expanded panel is meant to affect),
  * a net AMR shift: fraction of the genome painted AMR by panel A minus by B
    (does the expanded panel assign MORE Amerindigenous ancestry?),
  * per-position discordance so you can map WHERE panels disagree (e.g. near
    SLC16A11), for a figure.

If a truth is available (simulation), pair each panel against truth instead to
say which panel's differing tracts are actually CORRECT (use 02_score for that).
"""

import argparse
import itertools
import numpy as np
import pandas as pd
import pysam

AMR = 2  # ancestry index for AMR (matches truth_utils.ANC)


def load_flare_matrix(path):
    """FLARE .anc.vcf.gz -> (positions[M], calls[M, 2N]) with AN1/AN2 per sample."""
    vcf = pysam.VariantFile(path)
    samples = list(vcf.header.samples)
    pos, rows = [], []
    for rec in vcf:
        pos.append(rec.pos)
        r = []
        for s in samples:
            c = rec.samples[s]
            r.append(c.get("AN1")); r.append(c.get("AN2"))
        rows.append(r)
    return np.array(pos), np.array(rows, dtype=float), samples


def align(a_pos, a_call, b_pos, b_call):
    """Restrict both panels to shared positions (FLARE marker sets can differ)."""
    common, ia, ib = np.intersect1d(a_pos, b_pos, return_indices=True)
    return common, a_call[ia], b_call[ib]


def pair_stats(name_a, a_pos, a_call, name_b, b_pos, b_call):
    pos, A, B = align(a_pos, a_call, b_pos, b_call)
    mask = ~(np.isnan(A) | np.isnan(B))
    disc = (A != B) & mask
    overall = disc.sum() / mask.sum()
    a_amr = (A == AMR) & mask
    b_amr = (B == AMR) & mask
    amr_disagree = ((a_amr & ~b_amr) | (b_amr & ~a_amr)).sum() / mask.sum()
    net_amr_shift = (a_amr.sum() - b_amr.sum()) / mask.sum()
    per_pos = pd.DataFrame({
        "pos": pos,
        "discordance": disc.sum(axis=1) / np.maximum(mask.sum(axis=1), 1),
        "amr_shift": (a_amr.sum(axis=1) - b_amr.sum(axis=1)) / np.maximum(mask.sum(axis=1), 1),
    })
    summ = {
        "panel_A": name_a, "panel_B": name_b,
        "overall_discordance": round(float(overall), 4),
        "AMR_disagreement": round(float(amr_disagree), 4),
        "net_AMR_shift_A_minus_B": round(float(net_amr_shift), 4),
    }
    return summ, per_pos


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", action="append", required=True,
                    help="label=path.anc.vcf.gz (give 2+; all pairs compared)")
    ap.add_argument("--out_prefix", default="tract_differences")
    args = ap.parse_args()

    loaded = {}
    for spec in args.panel:
        label, path = spec.split("=", 1)
        p, c, _ = load_flare_matrix(path)
        loaded[label] = (p, c)

    summaries = []
    for a, b in itertools.combinations(loaded, 2):
        summ, per_pos = pair_stats(a, *loaded[a], b, *loaded[b])
        summaries.append(summ)
        per_pos.to_csv(f"{args.out_prefix}.{a}_vs_{b}.perpos.tsv.gz",
                       sep="\t", index=False, compression="gzip")

    df = pd.DataFrame(summaries)
    df.to_csv(f"{args.out_prefix}.summary.tsv", sep="\t", index=False)
    pd.set_option("display.width", 200)
    print(df.to_string(index=False))
    print("\nA positive net_AMR_shift (proposal minus baseline) means the expanded "
          "panel assigns MORE Amerindigenous ancestry; the per-pos files map where "
          "(overlay SLC16A11/ABCA1 for a figure). Pair vs simulated truth "
          "(02_score) to say which differing tracts are correct.")


if __name__ == "__main__":
    main()
