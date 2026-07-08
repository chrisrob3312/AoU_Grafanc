"""
Tier 2 (step 3) — simulate a phenotype with a KNOWN ancestry-specific effect,
using the SIMULATED individuals' TRUE local ancestry.

Because we control the truth, any difference in the recovered GWAS signal
between the two panels is attributable purely to LAI quality.

Model (per individual i, causal variant v):
    y_i = sum_a beta_a * (ALT copies of v on ancestry-a haplotypes) + noise
so beta_a is the per-ancestry effect. The headline setup gives AMR a real effect
and EUR/AFR ~0, mimicking an Amerindigenous-enriched signal (e.g. SLC16A11).

Outputs a phenotype file (person_id, y) + the ground-truth betas, for the
panel-vs-panel Tractor comparison (steps 4-5).

>>> This is the CORE experiment. <<<
"""

import argparse
import numpy as np
import pandas as pd
import pysam

from truth_utils import read_truth_bp, truth_at_positions


def alt_copies_by_ancestry(sim_vcf, truth_bp, causal_pos):
    """
    For each individual, count ALT copies of the causal variant on each ancestry
    background, using the phased sim genotypes + the true local ancestry at the
    causal position. Returns DataFrame[person_id, n_EUR, n_AFR, n_AMR].
    """
    truth = read_truth_bp(truth_bp)

    vcf = pysam.VariantFile(sim_vcf)
    samples = list(vcf.header.samples)
    rec = next(r for r in vcf if r.pos == causal_pos)

    counts = {s: {0: 0, 1: 0, 2: 0} for s in samples}
    for s in samples:
        alleles = rec.samples[s].allele_indices  # (a0, a1), phased
        for hap in (0, 1):
            if alleles[hap] != 1:
                continue  # REF on this haplotype, contributes nothing
            seg = truth.get((s, hap))
            if not seg:
                continue
            anc = truth_at_positions(seg, np.array([causal_pos]))[0]
            if anc in counts[s]:
                counts[s][anc] += 1
    return pd.DataFrame(
        [{"person_id": s, "n_EUR": c[0], "n_AFR": c[1], "n_AMR": c[2]}
         for s, c in counts.items()]
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sim_vcf", required=True, help="phased simulated genotypes")
    ap.add_argument("--truth_bp", required=True, help="haptools ground-truth breakpoints")
    ap.add_argument("--causal_pos", type=int, required=True, help="b38 pos of the causal variant")
    ap.add_argument("--beta_amr", type=float, default=0.5, help="AMR-specific effect")
    ap.add_argument("--beta_eur", type=float, default=0.0)
    ap.add_argument("--beta_afr", type=float, default=0.0)
    ap.add_argument("--noise_sd", type=float, default=1.0)
    ap.add_argument("--seed", type=int, default=1)   # explicit seed = reproducible
    ap.add_argument("--out_prefix", default="sim_pheno")
    args = ap.parse_args()

    counts = alt_copies_by_ancestry(args.sim_vcf, args.truth_bp, args.causal_pos)

    rng = np.random.default_rng(args.seed)
    genetic = (args.beta_eur * counts["n_EUR"]
               + args.beta_afr * counts["n_AFR"]
               + args.beta_amr * counts["n_AMR"])
    counts["y"] = genetic + rng.normal(0, args.noise_sd, len(counts))

    counts[["person_id", "y"]].to_csv(f"{args.out_prefix}.pheno.tsv", sep="\t", index=False)
    pd.DataFrame([{
        "causal_pos": args.causal_pos,
        "beta_EUR": args.beta_eur, "beta_AFR": args.beta_afr, "beta_AMR": args.beta_amr,
        "noise_sd": args.noise_sd, "seed": args.seed,
    }]).to_csv(f"{args.out_prefix}.truth.tsv", sep="\t", index=False)
    print(f"Wrote {args.out_prefix}.pheno.tsv and .truth.tsv  (n={len(counts)})")
    print(f"Injected per-ancestry betas: EUR={args.beta_eur} AFR={args.beta_afr} AMR={args.beta_amr}")


if __name__ == "__main__":
    main()
