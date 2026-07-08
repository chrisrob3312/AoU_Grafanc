# Evaluation — is the MXB-augmented AMR reference panel actually better?

Goal: show that an **expanded, well-matched Amerindigenous** reference panel
(1KG-HGDP + MX Biobank) gives better local-ancestry tracts for a **Latin
American** cohort — and better downstream discovery — than the alternatives,
**including the 1KG-style panel AoU will likely ship for CDR v9**.

## Design

- **Why FLARE for AoU:** it is far less cumbersome to run inside the AoU
  workbench than RFMix — a practical choice, not an accuracy/scaling claim.
- **Why RFMix1 appears at all:** the pre-AoU simulation accuracy benchmark uses
  **RFMix v1** on **admix-simu** truth so it is directly comparable to
  **Honorato-Mauer et al. 2024** (AJHG; bioRxiv 2024.08.26.609770 / PMC11866949).
  That study is the motivation: **AMR tracts are the weak spot** (TPR ~88–94% vs
  96–99% EUR/AFR), miscalls biased **AMR→EUR**, attributed to **small AMR
  reference size** — the gap MX Biobank fills.

### The core comparison: panels within FLARE

Do we get **different (and better) tracts** for a Latin American cohort with the
expanded MXB panel than with the alternatives? The arms (`EVAL_PANELS`)
deliberately **disentangle size vs matching**:

| Arm | AMR reference | Isolates |
|-----|---------------|----------|
| `mxb_expanded` | 1KG-HGDP + MXB, homogeneous | the proposal (large **and** matched) |
| `amr_small_homog` | few but pure Amerindigenous | matching alone (low N) |
| `amr_large_admixed` | many but admixed (1KG MXL/PEL/CLM) | size alone (contaminated) |
| `aou_default_1kg` | 1KG-style baseline | **what AoU CDR v9 will likely ship** |

`aou_default_1kg` is the key reference point — every result is framed as *what a
user gains over the tracts AoU itself provides*. And because Honorato-Mauer
blamed AMR **size**, the small-homog vs large-admixed vs MXB contrast shows
whether **size, matching, or both** is what actually matters.

Guiding principle: **real AoU data has no ground truth for local ancestry**, so
simulation gives the accuracy verdict, and real AoU data shows the
**consequences** — different tracts (step 8), stronger known-locus signal
(step 6), and candidate novel discovery (step 9). Test where the signal lives on
the **Amerindigenous** background, or a better AMR panel has nothing to show.

## Steps

| # | Question | Ground truth? | Files |
|---|----------|---------------|-------|
| 1 | Accuracy vs admix-simu truth — RFMix1 (Honorato-Mauer metrics) + FLARE, all panel arms; per-ancestry TPR and the AMR→EUR miscall asymmetry | Yes (admix-simu) | `01_simulate_admixed_truth.sh`, `02a_run_rfmix1.sh`, `02_score_lai_accuracy.py` |
| 8 | **Do the tracts actually differ across panels** on real AoU data (and where)? | No | `08_tract_differences.py` |
| 2 | Does better LAI improve a **Tractor GWAS** with a known injected ancestry-specific effect? | Yes (injected) | `03_..._phenotype.py`, `04_run_tractor_panelcompare.sh`, `05_gwas_power_bias.py` |
| 6 | **Known LAI regional association** — stronger AMR-track signal at real AMR-enriched loci? | No (published effects) | `06_realdata_positive_controls.py` |
| 9 | **Novel discovery** — AMR-track hits MXB finds that AoU-default misses (calibration-gated) | No | `09_novel_discovery.py` |
| 7 | Genome-wide calibration / admixture mapping / global-vs-local consistency | mixed | `07_robustness_and_value.py` |

**Story arc:** simulation accuracy (1) → tracts really change (8) → that changes
GWAS results causally (2) → it strengthens a known real signal (6) → and may
surface novel Latin American hits over the AoU-default panel (9), with
calibration guardrails throughout.

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
