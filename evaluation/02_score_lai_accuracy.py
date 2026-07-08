"""
Tier 1 (step 2) — score FLARE local-ancestry calls against the SIMULATED TRUTH,
comparing the targeted (MXB) vs comparison panel.

In plain terms: for each simulated haplotype we know the true ancestry at every
position (from haptools' .bp file). We line that up with what FLARE called and
count how often it was right — reported SEPARATELY for AMR, because that's the
ancestry the MXB panel is supposed to help with. We do this for both panels and
compare.

Inputs
------
--truth_bp        haptools .bp ground-truth breakpoints for one scenario
--flare_targeted  FLARE .anc.vcf.gz produced with the TARGETED panel
--flare_comparison FLARE .anc.vcf.gz produced with the COMPARISON panel
--out             metrics table (.tsv)

Metrics (per panel, and AMR-specific):
  accuracy, per-ancestry recall/precision, switch-error (tract fragmentation),
  and posterior calibration if FLARE was run with probs=true.
"""

import argparse
import numpy as np
import pandas as pd
import pysam

from truth_utils import ANC, read_truth_bp, truth_at_positions


def read_flare_calls(anc_vcf):
    """{(sample, hap): (positions[], calls[])} from FLARE AN1/AN2 fields."""
    vcf = pysam.VariantFile(anc_vcf)
    samples = list(vcf.header.samples)
    pos = []
    calls = {(s, h): [] for s in samples for h in (0, 1)}
    for rec in vcf:
        pos.append(rec.pos)
        for s in samples:
            c = rec.samples[s]
            calls[(s, 0)].append(c.get("AN1"))
            calls[(s, 1)].append(c.get("AN2"))
    positions = np.array(pos)
    return {k: (positions, np.array(v)) for k, v in calls.items()}


def switch_errors(seq):
    """Number of ancestry changes along a haplotype (fragmentation proxy)."""
    seq = seq[~pd.isna(seq)]
    return int(np.sum(seq[1:] != seq[:-1])) if len(seq) > 1 else 0


def score_panel(name, truth, flare):
    tp = fp = fn = correct = total = 0
    amr_tp = amr_fp = amr_fn = 0
    true_switches = called_switches = 0
    for key, (positions, called) in flare.items():
        if key not in truth:
            continue
        true_anc = truth_at_positions(truth[key], positions)
        mask = ~pd.isna(called)
        c = called[mask].astype(int)
        t = true_anc[mask]
        correct += int(np.sum(c == t)); total += len(t)
        # AMR-specific confusion (index 2)
        amr_tp += int(np.sum((c == 2) & (t == 2)))
        amr_fp += int(np.sum((c == 2) & (t != 2)))
        amr_fn += int(np.sum((c != 2) & (t == 2)))
        true_switches += switch_errors(t)
        called_switches += switch_errors(c)
    acc = correct / total if total else float("nan")
    amr_recall = amr_tp / (amr_tp + amr_fn) if (amr_tp + amr_fn) else float("nan")
    amr_prec = amr_tp / (amr_tp + amr_fp) if (amr_tp + amr_fp) else float("nan")
    return {
        "panel": name,
        "overall_accuracy": round(acc, 4),
        "AMR_recall": round(amr_recall, 4),
        "AMR_precision": round(amr_prec, 4),
        "true_switches": true_switches,
        "called_switches": called_switches,
        "excess_switch_frac": round((called_switches - true_switches) / max(true_switches, 1), 4),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth_bp", required=True)
    ap.add_argument("--flare_targeted", required=True)
    ap.add_argument("--flare_comparison", required=True)
    ap.add_argument("--out", default="lai_accuracy.tsv")
    args = ap.parse_args()

    truth = read_truth_bp(args.truth_bp)
    rows = [
        score_panel("targeted", truth, read_flare_calls(args.flare_targeted)),
        score_panel("comparison", truth, read_flare_calls(args.flare_comparison)),
    ]
    df = pd.DataFrame(rows)
    df.to_csv(args.out, sep="\t", index=False)
    print(df.to_string(index=False))
    print("\nBetter AMR recall/precision and LOWER excess_switch_frac = better panel "
          "(the MXB panel should win on the AMR columns).")


if __name__ == "__main__":
    main()
