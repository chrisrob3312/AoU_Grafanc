#!/usr/bin/env bash
# Configuration for the reference-panel evaluation.
#   source config/config.sh; source evaluation/eval_config.sh
set -euo pipefail

: "${GA_ROOT:?source config/config.sh first}"

export EVAL_ROOT="${GA_ROOT}/evaluation"
export EVAL_SIM_DIR="${EVAL_ROOT}/sim"          # simulated cohorts + ground truth
export EVAL_LAI_DIR="${EVAL_ROOT}/lai"          # FLARE calls on simulated data
export EVAL_GWAS_DIR="${EVAL_ROOT}/gwas"        # Tractor runs + metrics
export EVAL_REAL_DIR="${EVAL_ROOT}/realdata"    # real positive-control loci

# ---- Simulation parameters (Tier 1/2) --------------------------------------
export SIM_N_INDIV=5000            # simulated admixed individuals per scenario
export SIM_GENERATIONS=10          # generations since admixture (tune to cohort)
export SIM_CHROM="chr17"           # chromosome to simulate (SLC16A11 is chr17)
# Admixture proportions per scenario "label:AFR,EUR,AMR"
export SIM_SCENARIOS=("afr_eur:0.5,0.5,0.0" "afr_eur_amr:0.25,0.45,0.30")

# ---- Two-track evaluation design -------------------------------------------
# TRACK A (accuracy yardstick): admix-simu + RFMix1, scored with the metrics of
#   Honorato-Mauer et al. 2024 (AJHG; "Characterizing features affecting local
#   ancestry inference performance in admixed populations") so our AMR-accuracy
#   numbers are directly comparable to that prior benchmark. That paper shows
#   AMR tracts are the weak spot (TPR ~88-94% vs 96-99% EUR/AFR) with miscalls
#   biased AMR->EUR, attributed to the SMALL AMR reference size -- exactly the
#   gap MX Biobank fills.
# TRACK B (biobank scale): FLARE for the actual AoU run (RFMix1 does not scale to
#   AoU N). A FLARE-vs-RFMix1 concordance on the simulated data shows FLARE
#   reproduces the RFMix1 accuracy while scaling.

# ---- Tools -----------------------------------------------------------------
# admix-simu (Williams/Martin): simulate admixed haplotypes WITH ground-truth
# breakpoints. Run LOCALLY before the AoU work; we ingest its outputs here.
export ADMIXSIMU_DIR="${ADMIXSIMU_DIR:-${HOME}/tools/admix-simu}"   # has simu-mix.pl
# RFMix v1 (Track A accuracy benchmark, matches Honorato-Mauer).
export RFMIX1_DIR="${RFMIX1_DIR:-${HOME}/tools/RFMix_v1.5.4}"       # has RunRFMix.py
export TRACTOR_DIR="${HOME}/tools/Tractor"      # git clone of Atkinson-Lab/Tractor

# ---- Positive-control loci (Tier 3). GRCh38. --------------------------------
# name | chrom | pos (b38) | trait | expected ancestry | role
export EVAL_LOCI=(
  "SLC16A11|chr17|7041768|T2D|AMR|amr_primary"
  "ABCA1_R230C|chr9|104903697|HDL|AMR|amr_secondary"
  "HNF1A_E508K|chr12|120994371|T2D|AMR|amr_secondary"
  "DARC_rs2814778|chr1|159205564|neutrophil|AFR|specificity_control"
  "SLC24A5|chr15|48134287|pigmentation|EUR_AFR|specificity_control"
  "APOL1|chr22|36265860|kidney|AFR|specificity_control"
)
# NOTE: confirm each b38 position against dbSNP for the exact tag SNP you use.

# ---- Real-data trait anchor (Tier 3) ---------------------------------------
# Default: Type 2 diabetes via SLC16A11.  Pull the phenotype from the CDR.
export REAL_TRAIT="T2D"
# OMOP concept ids / SQL for the case definition (fill from your CDR).
export REAL_TRAIT_SQL="<BigQuery SQL selecting person_id, phenotype from ${WORKSPACE_CDR}>"

echo "Loaded eval config → ${EVAL_ROOT}"
