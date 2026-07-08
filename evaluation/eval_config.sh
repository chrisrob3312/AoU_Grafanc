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

# ---- Evaluation design ------------------------------------------------------
# WHY FLARE for AoU: it is far less cumbersome to run inside the AoU workbench
# than RFMix -- that is the practical reason, NOT a scaling/accuracy claim.
#
# WHY RFMix1 at all: the pre-AoU simulation accuracy benchmark uses RFMix v1 so
# it is directly comparable to Honorato-Mauer et al. 2024 (AJHG; PMC11866949),
# which showed AMR/Native American tracts are the weak spot (TPR ~88-94% vs
# 96-99% EUR/AFR), miscalls biased AMR->EUR, attributed to SMALL AMR reference
# size. That is the gap we test.
#
# THE CORE QUESTION: within FLARE, do we get DIFFERENT (and better) tracts for a
# Latin American cohort with the expanded MXB panel than with the alternatives,
# and does that translate into downstream discovery power? The alternatives are
# defined in EVAL_PANELS below and deliberately disentangle SIZE vs MATCHING:
#   * a small but homogeneous Amerindigenous AMR reference (matched, low N),
#   * a larger but admixed / less-matched AMR reference (high N, contaminated),
#   * the expanded MXB reference (high N AND homogeneous -- the proposal),
#   * an AoU-default-style 1KG reference: AoU CDR v9 will release LAI tracts, but
#     likely WITHOUT an expanded Amerindigenous panel (probably 1KG-based) -- so
#     this arm shows what a user gains over the tracts AoU itself will ship.
#
# label | reference-sample keep-list (built in flare/ref_panel step 8 variants)
export EVAL_PANELS=(
  "mxb_expanded|<gs://.../ref/mxb_expanded.keep>"       # 1KG-HGDP + MXB, homogeneous AMR (proposal)
  "amr_small_homog|<gs://.../ref/amr_small_homog.keep>" # few but pure Amerindigenous AMR
  "amr_large_admixed|<gs://.../ref/amr_large_admixed.keep>" # many but admixed AMR (e.g. 1KG MXL/PEL/CLM)
  "aou_default_1kg|<gs://.../ref/aou_default_1kg.keep>" # 1KG-style baseline ~ what AoU v9 will ship
)
# Which arm is the reference point for "improvement over what AoU ships":
export EVAL_BASELINE_PANEL="aou_default_1kg"
# Which arm is the proposal:
export EVAL_PROPOSAL_PANEL="mxb_expanded"

# ---- Tools -----------------------------------------------------------------
# admix-simu (Williams/Martin): simulate admixed haplotypes WITH ground-truth
# breakpoints. Run LOCALLY before the AoU work; we ingest its outputs here.
export ADMIXSIMU_DIR="${ADMIXSIMU_DIR:-${HOME}/tools/admix-simu}"   # has simu-mix.pl
# RFMix v1 (simulation accuracy benchmark, matches Honorato-Mauer).
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
