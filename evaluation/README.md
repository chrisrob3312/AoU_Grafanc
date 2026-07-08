# Evaluation — is the MXB-augmented AMR reference panel actually better?

Goal: quantify the value of the **targeted** (1KG-HGDP + MX Biobank, >95%
homogeneous AMR) reference panel over a **comparison** (broad/admixed AMR, no
MXB) panel for local-ancestry inference (LAI) and the downstream Tractor GWAS.

## Two-track design

- **Track A — accuracy yardstick: admix-simu + RFMix1.** We benchmark LAI
  accuracy with **RFMix v1** on **admix-simu** ground truth, using the metrics of
  **Honorato-Mauer et al. 2024** (AJHG, *Characterizing features affecting local
  ancestry inference performance in admixed populations*; bioRxiv
  2024.08.26.609770 / PMC11866949). That study is the motivation: **AMR/Native
  American tracts are the weak spot** (TPR ~88–94% vs 96–99% EUR/AFR) with
  miscalls biased **AMR→EUR**, which they attribute to the **small AMR reference
  sample size** — exactly the gap MX Biobank fills. So the headline Track-A
  result is: does adding MXB **raise AMR TPR and shrink the AMR→EUR miscall
  asymmetry**, comparably measured to that paper?
- **Track B — biobank scale: FLARE.** RFMix1 does not scale to AoU N; FLARE does.
  We run FLARE on the same simulated data and on real AoU, and check
  **FLARE-vs-RFMix1 concordance** on the sims to show FLARE reproduces the
  RFMix1 accuracy while scaling. The improved MXB panel is what FLARE then
  carries at biobank scale.

Guiding principle: **real AoU data has no ground truth for local ancestry**, so
the cleanest evidence is simulation; real positive-control loci are
corroboration. And you must test at loci where the signal lives on the
**Amerindigenous** background, or a better AMR panel has nothing to show.

## Four tiers

| Tier | Question | Ground truth? | Files |
|------|----------|---------------|-------|
| 1 | Does the MXB panel call AMR **tracts** more accurately (RFMix1 vs Honorato-Mauer metrics; FLARE-vs-RFMix1 concordance)? | Yes (admix-simu) | `01_simulate_admixed_truth.sh`, `02a_run_rfmix1.sh`, `02_score_lai_accuracy.py` |
| 2 | Does better LAI improve a **Tractor GWAS** with a known ancestry-specific effect? | Yes (injected) | `03_simulate_ancestry_specific_phenotype.py`, `04_run_tractor_panelcompare.sh`, `05_gwas_power_bias.py` |
| 3 | Does it reproduce on **real** AMR-enriched loci? | No (published effects) | `06_realdata_positive_controls.py` |
| 4 | Genome-wide **calibration / robustness / extra value** | mixed | `07_robustness_and_value.py` |

**Lead result** = Tier 2 (panel-vs-panel Tractor on a simulated ancestry-specific
phenotype), backed by Tier 1 LAI accuracy, validated on 2–3 Tier-3 loci.

## Positive-control loci (config in `eval_config.sh`)

| Locus | Trait | Ancestry-specific signal | Role |
|-------|-------|--------------------------|------|
| **SLC16A11** | Type 2 diabetes | Risk haplotype common in Native American/Latino, ~absent AFR, rare EUR (SIGMA, *Nature* 2014) | **AMR primary** |
| **ABCA1** rs9282541 (R230C) | Low HDL | Essentially private to Amerindigenous/Mestizo | AMR secondary |
| **HNF1A** E508K | T2D (Latino) | Latino-enriched | AMR secondary |
| **DARC/ACKR1** rs2814778 | Neutrophil count | Near-fixed AFR vs others | Specificity / negative-for-AMR control |
| **SLC24A5** | Pigmentation | EUR/AFR | Specificity control |
| **APOL1** G1/G2 | Kidney | AFR | Specificity control |

The AMR loci should gain power / de-bias with the MXB panel; the AFR/EUR
controls should NOT change with the AMR panel (that's the specificity argument).

## Critical design guard

The **simulation founders must be held out of the FLARE reference panel** — if
the same individuals are used to both simulate and paint, accuracy is inflated.
`01_simulate_admixed_truth.sh` splits founders vs panel by sample ID.

## Why Tractor makes the AMR panel visible

Tractor deconvolves association by local ancestry, estimating a **per-ancestry
effect and allele count**. If AMR tracts are mis-assigned, AMR alleles leak into
the EUR/AFR bins → the AMR-specific effect **attenuates and mislocates**. A
better-matched AMR panel sharpens and de-biases that AMR estimate — which is
exactly what Tiers 2–3 measure.
