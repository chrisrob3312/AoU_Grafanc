# AoU v9 GrafAnc → FLARE local-ancestry pipeline

Global ancestry (GrafAnc) for **every** All of Us CDR **v9** participant with
imputed-array **or** WGS data, then FLARE local-ancestry inference (LAI) for the
admixed cohorts, post-processed for [Tractor](https://github.com/Atkinson-Lab/Tractor).
Built for the Verily Researcher Workbench.

## What it does

1. **GrafAnc global ancestry** on the union of imputed + WGS participants, using
   the 282,424 GrafAnc ancestry SNPs, with a light in-house QC pass (no
   lab-curated QC list needed for v9).
2. **Categorize** by `AncGroupID`: African American (107), Latin American 1/2
   (601/602 = "LA1"/"LA2"), Native American (603), Multiracial (800).
3. **Resolve Multiracial** with a supervised 7-superpop ADMIXTURE run
   (AMR/AFR/EUR/SAS/EAS/MEN/OCN) to pull out true 3-way AFR-EUR-AMR individuals.
4. **FLARE LAI**, 2-way (AFR-EUR) and 3-way (AFR-EUR-AMR), against two reference
   panels — your **targeted** 1KG-HGDP + MXB panel and a **comparison** panel —
   to benchmark panel fit.
5. **Post-process**: propagate calls to uncalled variants by nearest tract, then
   format for Tractor.

## Run order

| Step | File | Env |
|------|------|-----|
| — | `config/config.sh` | edit first; sourced by every script |
| 1 | `scripts/01_setup_grafanc.sh` | RW terminal |
| 2 | `notebooks/02_extract_anc_snps.py` | Hail / Dataproc |
| 3 | `notebooks/03_qc_and_export.py` | Hail / Dataproc |
| 4 | `scripts/04_run_grafanc.sh` | RW terminal |
| 5 | `notebooks/05_categorize_cohorts.py` | Jupyter (pandas) |
| 6 | `admixture/06_multiracial_supervised_admixture.sh` | RW terminal (PLINK + ADMIXTURE) |
| 7 | `notebooks/07_define_flare_cohorts.py` | Jupyter (pandas) |
| 8 | `flare/ref_panel/08_build_ref_panels.sh` | RW terminal (bcftools + SHAPEIT5) |
| 9 | `flare/09_phase_target.sh <cohort> <prefix>` | RW terminal |
| 10 | `flare/10_run_flare.sh <2way\|3way> <targeted\|comparison>` | RW terminal (Java) |
| 11 | `flare/postprocess/11_fill_uncalled_tracts.py` | Jupyter |
| 12 | `flare/postprocess/12_flare_to_tractor.py` | Jupyter |
| 13 | `flare/compare/13_panel_accuracy_comparison.py` | Jupyter |

## Cohort → FLARE mode map

| GrafAnc group | AncGroupID | FLARE |
|---|---|---|
| African American | 107 | 2-way AFR-EUR |
| Latin American 1 / 2 | 601 / 602 | 3-way AFR-EUR-AMR |
| Native American | 603 | AMR anchor (ref / eval) |
| Multiracial (AFR-EUR-AMR) | 800 → step 6/7 | 3-way AFR-EUR-AMR |
| Multiracial (near AFR-EUR) | 800 → step 6/7 | 2-way AFR-EUR |

## What you still need to fill in

Search the tree for `<...>` placeholders. The load-bearing ones:

- **`config/config.sh`** — workspace name, v9 genomic env-var names/paths
  (WGS VDS, imputed MT, ACAF callset), Dataproc worker count, FLARE heap/mem.
- **ADMIXTURE homogeneous sample lists** (`REF_HOMOG_{EUR,AFR,AMR}_SAMPLES`,
  `REF_COMPARISON_AMR_SAMPLES`) — from your local HPC ADMIXTURE run. *Pending.*
- **Combined reference callset** path (`REF_COMBINED_VCF`) and the 7-superpop
  reference/`.pop` for step 6.
- **Genetic maps** (`GENETIC_MAP_DIR`) for phasing + FLARE.
- **Lab templates** — `flare_to_tracts()` / `assign_variant_ancestry()` in
  step 11 and `write_tractor_matrices()` in step 12 have runnable defaults;
  swap in the lab versions when they arrive. *Pending.*
- **FLARE FORMAT field names** (`AN1`/`AN2`/`ANP`) — confirm against the FLARE
  version you download; they're flagged inline.

See `docs/PIPELINE.md` for the design rationale, QC choices, and open decisions.

## Tools

- GrafAnc — <https://github.com/jimmy-penn/grafanc> ([paper](https://www.sciencedirect.com/science/article/pii/S2666247725001332))
- FLARE — <https://github.com/browning-lab/flare>
- Tractor — <https://github.com/Atkinson-Lab/Tractor>
