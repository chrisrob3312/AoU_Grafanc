"""
Tier 1 (step 2) — score LAI calls against the admix-simu ground truth, for BOTH
methods (RFMix1 = accuracy yardstick; FLARE = biobank-scale) and BOTH panels
(targeted MXB vs comparison), in one table.

Metrics follow Honorato-Mauer et al. 2024 so our AMR numbers are comparable:
  * per-ancestry TPR (true-positive rate / recall), reported for EUR/AFR/AMR;
  * the AMR<->EUR MISCALL ASYMMETRY (fraction of true-AMR called EUR vs true-EUR
    called AMR) — the paper's key finding, driven by small AMR reference size,
    which MX Biobank is meant to fix.

Pass any number of call sets:
  --flare   label=path.anc.vcf.gz            (repeatable)
  --rfmix1  label=viterbi.txt:positions.txt:ids.txt   (repeatable)
  --rfmix_classes EUR,AFR,AMR                (1-based class order in the .classes)

Rows for {rfmix1,flare} x {targeted,comparison} let you read off (a) MXB gain in
AMR TPR within each method and (b) FLARE-vs-RFMix1 concordance for scaling.
"""

import argparse
import numpy as np
import pandas as pd
import pysam

from truth_utils import ANC, read_truth_bp, truth_at_positions


def load_flare(path):
    """FLARE .anc.vcf.gz -> {(sample,hap):(positions[], calls[])} via AN1/AN2."""
    vcf = pysam.VariantFile(path)
    samples = list(vcf.header.samples)
    pos, calls = [], {(s, h): [] for s in samples for h in (0, 1)}
    for rec in vcf:
        pos.append(rec.pos)
        for s in samples:
            c = rec.samples[s]
            calls[(s, 0)].append(c.get("AN1"))
            calls[(s, 1)].append(c.get("AN2"))
    positions = np.array(pos)
    return {k: (positions, np.array(v, dtype=float)) for k, v in calls.items()}


def load_rfmix1(viterbi, positions_file, ids_file, class_order):
    """
    RFMix1 Viterbi -> {(sample,hap):(positions[], calls[])}.
    Viterbi: rows = markers/windows, cols = haplotypes (2 per individual, in
    query order), values = 1-based ancestry class. class_order maps class i to a
    label; we convert to the ANC index used by the truth.
    """
    vit = np.loadtxt(viterbi, dtype=int)                 # markers x haplotypes
    positions = np.loadtxt(positions_file, dtype=int)    # bp per marker row
    ids = [l.strip() for l in open(ids_file) if l.strip()]
    label_of_class = {i + 1: lab for i, lab in enumerate(class_order)}
    anc_idx = {lab: idx for idx, lab in ANC.items()}
    out = {}
    for j in range(vit.shape[1]):
        sample = ids[j // 2]
        hap = j % 2
        mapped = np.array([anc_idx.get(label_of_class.get(v, ""), np.nan) for v in vit[:, j]])
        out[(sample, hap)] = (positions, mapped)
    return out


def score(label, method, truth, calls):
    # confusion counts over all aligned (position, haplotype) calls
    conf = np.zeros((3, 3), dtype=np.int64)   # conf[true, called]
    for key, (positions, called) in calls.items():
        if key not in truth:
            continue
        t = truth_at_positions(truth[key], positions)
        mask = ~pd.isna(called)
        c = called[mask].astype(int)
        tt = t[mask]
        for a in range(3):
            for b in range(3):
                conf[a, b] += int(np.sum((tt == a) & (c == b)))
    tpr = {ANC[a]: (conf[a, a] / conf[a].sum() if conf[a].sum() else float("nan"))
           for a in range(3)}
    # AMR(2)<->EUR(0) miscall asymmetry (paper's headline)
    amr_to_eur = conf[2, 0] / conf[2].sum() if conf[2].sum() else float("nan")
    eur_to_amr = conf[0, 2] / conf[0].sum() if conf[0].sum() else float("nan")
    return {
        "label": label, "method": method,
        "TPR_EUR": round(tpr["EUR"], 4), "TPR_AFR": round(tpr["AFR"], 4),
        "TPR_AMR": round(tpr["AMR"], 4),
        "AMR->EUR_miscall": round(amr_to_eur, 4),
        "EUR->AMR_miscall": round(eur_to_amr, 4),
        "AMR_EUR_asymmetry": round(amr_to_eur - eur_to_amr, 4),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth_bp", required=True)
    ap.add_argument("--truth_ids", required=True)
    ap.add_argument("--flare", action="append", default=[], help="label=path.anc.vcf.gz")
    ap.add_argument("--rfmix1", action="append", default=[],
                    help="label=viterbi.txt:positions.txt:ids.txt")
    ap.add_argument("--rfmix_classes", default="EUR,AFR,AMR")
    ap.add_argument("--out", default="lai_accuracy.tsv")
    args = ap.parse_args()

    truth = read_truth_bp(args.truth_bp, args.truth_ids)
    class_order = args.rfmix_classes.split(",")
    rows = []
    for spec in args.flare:
        label, path = spec.split("=", 1)
        rows.append(score(label, "flare", truth, load_flare(path)))
    for spec in args.rfmix1:
        label, triple = spec.split("=", 1)
        vit, pos, ids = triple.split(":")
        rows.append(score(label, "rfmix1", truth, load_rfmix1(vit, pos, ids, class_order)))

    df = pd.DataFrame(rows)
    df.to_csv(args.out, sep="\t", index=False)
    pd.set_option("display.width", 200)
    print(df.to_string(index=False))
    print("\nReading it:")
    print("  * Within a method, targeted (MXB) should raise TPR_AMR and SHRINK "
          "AMR_EUR_asymmetry vs comparison — the Honorato-Mauer gap closing.")
    print("  * rfmix1 vs flare rows on the same panel = concordance justifying "
          "FLARE for biobank-scale AoU.")


if __name__ == "__main__":
    main()
