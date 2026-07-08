#!/usr/bin/env bash
# Tier 1 (step 2a, Track A) — run RFMix v1 on the admix-simu cohorts.
#
# RFMix1 is our accuracy yardstick because it matches Honorato-Mauer et al. 2024.
# Run it with the TARGETED (MXB) panel and the COMPARISON panel so we can measure
# whether adding MXB raises AMR tract accuracy and reduces the AMR->EUR miscall
# asymmetry that paper attributes to small AMR reference size.
#
# RFMix1 does NOT scale to AoU N — it is used ONLY on the simulated benchmark.
# FLARE (Track B) carries the same panels at AoU biobank scale.
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"
source "$(dirname "$0")/eval_config.sh"

SCENARIO="${1:?usage: 02a_run_rfmix1.sh <scenario_label>}"
WORK="${HOME}/eval_rfmix1/${SCENARIO}"; mkdir -p "${WORK}"; cd "${WORK}"

# Held-out reference (founders removed in step 1) + simulated targets.
gsutil -m cp "${EVAL_SIM_DIR}/eval_panel.heldout.vcf.gz"* .
gsutil -m cp "${EVAL_SIM_DIR}/${SCENARIO}/"* .

for PANEL in targeted comparison; do
  # RFMix1 classes come from the panel's ancestry labels (EUR/AFR/AMR); the
  # comparison panel swaps in the broad-AMR (no-MXB) reference samples.
  REF_KEEP="${FLARE_REF_DIR}/${PANEL}/${PANEL}.keep"     # reused from panel build
  gsutil cp "${REF_KEEP}" ./${PANEL}.keep

  # RFMix1 wants: phased alleles (0/1), a classes file (ref pop per hap; 0 for
  # query/admixed), and a genetic-position (cM) file. Build them with the helper
  # scripts shipped in RFMix1 / ancestry_pipeline, or bcftools + a small awk.
  python "${RFMIX1_DIR}/scripts/make_rfmix_inputs.py" \
      --query   "sim_${SCENARIO}.phased.vcf.gz" \
      --ref     "eval_panel.heldout.vcf.gz" \
      --ref_keep "${PANEL}.keep" \
      --ref_labels "<gs://.../7superpop_ref_labels.tsv>" \
      --genetic_map "<${GENETIC_MAP_DIR}/${SIM_CHROM}.cM.map>" \
      --out_prefix "rfmix_in_${PANEL}"

  python "${RFMIX1_DIR}/RunRFMix.py" PopPhased \
      "rfmix_in_${PANEL}.alleles" \
      "rfmix_in_${PANEL}.classes" \
      "rfmix_in_${PANEL}.snp_locations" \
      -o "rfmix_out_${PANEL}" \
      --num-threads "${N_THREADS}"
  # rfmix_out_${PANEL}.*.Viterbi.txt = per-window ML ancestry per haplotype.

  gsutil -m cp "rfmix_out_${PANEL}".* "rfmix_in_${PANEL}.snp_locations" \
      "${EVAL_LAI_DIR}/${SCENARIO}/"
done
echo "RFMix1 done for ${SCENARIO}. Score with 02_score_lai_accuracy.py --method rfmix1"
