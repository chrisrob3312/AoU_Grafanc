#!/usr/bin/env bash
# Tier 1 (step 1) — obtain admix-simu simulated individuals WITH ground-truth
# local ancestry, and stage them for LAI.
#
# You run admix-simu LOCALLY before the AoU work (that's your existing pipeline,
# matched to Honorato-Mauer et al. 2024). This script (a) ingests those outputs
# into the workspace bucket, and (b) enforces the critical guard: the simulation
# FOUNDERS must be held out of the LAI reference panel, or accuracy is inflated.
# It also documents the admix-simu invocation for reproducibility.
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"
source "$(dirname "$0")/eval_config.sh"

WORK="${HOME}/eval_sim"; mkdir -p "${WORK}"; cd "${WORK}"

# ---- Inputs from your local admix-simu run ----------------------------------
# For each scenario: the simulated phased haplotypes, the .bp ground truth, and
# the sample-order (.ids) file that names the individuals in .bp order.
LOCAL_SIM_PREFIX="<gs://.../admixsimu/{label}/sim_{label}>"   # .hap/.vcf.gz, .bp, .ids
# The founder sample IDs admix-simu drew from (so we can hold them out).
SIM_FOUNDERS="<gs://.../admixsimu/founders.samples.txt>"

# Reference build for admix-simu (documented for reproducibility; run locally):
#   ${ADMIXSIMU_DIR}/simu-mix.pl model_{label}.dat map.txt sim_{label} \
#       --founder_haps founders.phgeno --founder_info founders.ids
# where model_{label}.dat encodes SIM_GENERATIONS + the per-pop admixture
# fractions in SIM_SCENARIOS.

gsutil cp "${SIM_FOUNDERS}" ./founders.keep

for scen in "${SIM_SCENARIOS[@]}"; do
  label="${scen%%:*}"
  pfx="${LOCAL_SIM_PREFIX//\{label\}/${label}}"
  gsutil -m cp "${pfx}".* "${EVAL_SIM_DIR}/${label}/"
  echo "Ingested admix-simu scenario '${label}' → ${EVAL_SIM_DIR}/${label}/"
done

# ---- Hold founders out of the LAI reference panel ---------------------------
# Build the panel-side reference = combined callset MINUS the simulation founders.
gsutil -m cp "${REF_COMBINED_VCF}"* .
COMBINED="$(basename "${REF_COMBINED_VCF}")"
bcftools view -S ^founders.keep --force-samples "${COMBINED}" \
  -Oz -o eval_panel.heldout.vcf.gz
bcftools index -t eval_panel.heldout.vcf.gz
gsutil -m cp eval_panel.heldout.vcf.gz* founders.keep "${EVAL_SIM_DIR}/"

echo "Held-out reference panel → ${EVAL_SIM_DIR}/eval_panel.heldout.vcf.gz"
echo "NEXT (Track A): RFMix1 on the sims  -> 02a_run_rfmix1.sh"
echo "     (Track B): FLARE on the sims   -> main flow steps 9-10 (restrict ref to held-out panel)"
echo "     then score both -> 02_score_lai_accuracy.py"
