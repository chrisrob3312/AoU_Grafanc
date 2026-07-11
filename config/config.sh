#!/usr/bin/env bash
# =============================================================================
# Central configuration for the AoU v9 GrafAnc + FLARE local-ancestry pipeline.
#
#   source config/config.sh
#
# Everything the pipeline needs is defined here so the individual scripts stay
# portable across workspaces.  Anything wrapped in <> you must fill in.
# AoU-provided environment variables are read from the workbench and should NOT
# be hard-coded — we only fall back to a documented default when unset.
# =============================================================================

set -euo pipefail

# ---- Workspace identity -----------------------------------------------------
export WORKSPACE_NAME="<copy and paste workspace name>"
export CDR_VERSION="v9"          # All of Us Controlled Tier CDR release

# ---- AoU-provided environment variables (Verily Researcher Workbench) -------
# These are exported into every RW analysis environment.  If your v9 workbench
# uses different names, confirm them under "Genomic data" in the workbench docs
# and override here.  The values below are the v8 names as a fallback.
: "${WORKSPACE_BUCKET:?WORKSPACE_BUCKET not set — are you inside the RW?}"
: "${WORKSPACE_CDR:?WORKSPACE_CDR (BigQuery dataset) not set}"

# WGS short-read joint callset (Hail VDS)
export WGS_VDS_PATH="${WGS_VDS_PATH:-<gs://.../wgs/.../vds/hail.vds>}"

# Imputed genotype callset.  In AoU this is the array-based imputed data.
# Confirm the v9 path/format (Hail MT, VDS, or per-chrom BGEN/VCF).
export IMPUTED_MT_PATH="${IMPUTED_MT_PATH:-<gs://.../imputation/.../hail.mt>}"

# ---- Where our outputs live (all under the workspace bucket) ----------------
export GA_ROOT="${WORKSPACE_BUCKET}/grafanc_flare_${CDR_VERSION}"
export GA_SNP_DIR="${GA_ROOT}/anc_snps"           # extracted ancestry SNPs
export GA_QC_DIR="${GA_ROOT}/qc"                  # post-QC genotypes for GrafAnc
export GA_RESULTS_DIR="${GA_ROOT}/grafanc_results"
export GA_COHORT_DIR="${GA_ROOT}/cohorts"         # per-ancestry sample lists
export FLARE_ROOT="${GA_ROOT}/flare"
export FLARE_REF_DIR="${FLARE_ROOT}/ref_panels"
export FLARE_OUT_DIR="${FLARE_ROOT}/output"
export FLARE_POST_DIR="${FLARE_ROOT}/postprocess"

# ---- Tool source repositories ----------------------------------------------
export GRAFANC_REPO="https://github.com/jimmy-penn/grafanc.git"
export FLARE_JAR_URL="https://faculty.washington.edu/browning/flare.jar"   # Browning lab
export TRACTOR_REPO="https://github.com/Atkinson-Lab/Tractor.git"

# GrafAnc ships its 282,424 ancestry SNPs in cpp/data/AncSnpPopAFs.txt
# (col1 chrom, col2 GRCh37 pos, col3 GRCh38 pos, col4 RS ID).
export GRAFANC_DIR="${HOME}/tools/grafanc"
export ANC_SNP_TABLE="${GRAFANC_DIR}/cpp/data/AncSnpPopAFs.txt"

# ---- GrafAnc ancestry-group IDs we care about -------------------------------
# (First digit = continental super-pop; see docs.)
#   107 = African American        601 = Latin American 1 ("LA1")
#   602 = Latin American 2 ("LA2")  603 = Native American
#   800 = Multiracial (resolved further by supervised ADMIXTURE)
export ANCGROUP_AFRICAN_AMERICAN=107
export ANCGROUP_LATIN_AMERICAN_1=601
export ANCGROUP_LATIN_AMERICAN_2=602
export ANCGROUP_NATIVE_AMERICAN=603
export ANCGROUP_MULTIRACIAL=800

# ---- QC thresholds for the ancestry SNPs ------------------------------------
# GrafAnc tolerates missingness and needs no LD pruning, so QC is deliberately
# light: quality + call rate + allele sanity only.  We do NOT apply a global
# HWE filter — the cohort is genetically structured, so HWE is violated by
# ancestry, not by genotyping error.
export QC_MIN_CALL_RATE=0.95      # per-variant across the whole cohort
export QC_WGS_MIN_GQ=20           # WGS genotype quality floor
export QC_WGS_MIN_DP=10           # WGS read depth floor
export QC_IMPUTED_MIN_R2=0.80     # imputation quality (INFO/R2) floor

# ---- FLARE reference-panel definitions --------------------------------------
# Combined joint-called 1KG-HGDP + 50 MX Biobank WGS callset (your local build).
export REF_COMBINED_VCF="<gs://.../1kg_hgdp_mxb.joint.vcf.gz>"

# Homogeneous per-super-pop reference lists from the local ADMIXTURE homogeneity
# screen — committed in reference/homog_by_super_pop/ (AFR 634, AMR 88 incl. 50
# MXB, EAS 667, EUR 620, SAS 48). Copy these into the workspace bucket, or point
# at the repo copies. AMR (88) IS the expanded/MXB panel; drop MXB_* for the
# small-homogeneous arm (reference/build_panels.sh does the split).
REF_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../reference" 2>/dev/null && pwd || echo '<repo>/reference')"
export REF_HOMOG_EUR_SAMPLES="${REF_HOMOG_EUR_SAMPLES:-${REF_DIR}/homog_by_super_pop/EUR_ids.txt}"
export REF_HOMOG_AFR_SAMPLES="${REF_HOMOG_AFR_SAMPLES:-${REF_DIR}/homog_by_super_pop/AFR_ids.txt}"
export REF_HOMOG_AMR_SAMPLES="${REF_HOMOG_AMR_SAMPLES:-${REF_DIR}/homog_by_super_pop/AMR_ids.txt}"

# Comparison ("less-fitting") AMR panel: the 38 non-MXB Amerindigenous from the
# homogeneity panel (matched but small). For the larger/admixed and 1KG-default
# arms, supply external AMR lists (see reference/README.md + evaluation/).
export REF_COMPARISON_AMR_SAMPLES="${REF_COMPARISON_AMR_SAMPLES:-${REF_DIR}/panels/AMR_nonmxb.txt}"

# Genetic map for FLARE / phasing (GRCh38, PLINK format cM map per chrom).
export GENETIC_MAP_DIR="<gs://.../genetic_maps/plink.GRCh38>"

# ---- Compute defaults -------------------------------------------------------
export N_THREADS="${N_THREADS:-8}"
export N_WORKERS="${N_WORKERS:-<num_dataproc_workers>}"

echo "Loaded config for workspace '${WORKSPACE_NAME}' (CDR ${CDR_VERSION})."
echo "  Outputs → ${GA_ROOT}"
