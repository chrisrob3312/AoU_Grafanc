#!/usr/bin/env bash
# Tier 2 (step 4) — run the Tractor GWAS on the simulated phenotype TWICE:
# once using LAI from the targeted panel, once from the comparison panel.
# Everything else (genotypes, phenotype, covariates) is held fixed, so any
# difference in the GWAS is attributable to LAI quality alone.
#
# Pipeline per panel:
#   FLARE .anc.vcf.gz  --(step 11)-->  per-variant ancestry
#                      --(step 12)-->  Tractor dosage/hapcount matrices
#   Tractor ancestry-aware regression on those matrices + the phenotype.
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"
source "$(dirname "$0")/eval_config.sh"

SCENARIO="${1:?usage: 04_run_tractor_panelcompare.sh <scenario_label> <pheno.tsv>}"
PHENO="${2:?pass the simulated phenotype file from step 3}"

[[ -d "${TRACTOR_DIR}" ]] || git clone "${TRACTOR_REPO}" "${TRACTOR_DIR}"

WORK="${HOME}/eval_gwas/${SCENARIO}"; mkdir -p "${WORK}"; cd "${WORK}"
gsutil cp "${PHENO}" ./pheno.tsv

for PANEL in targeted comparison; do
  echo "=== Tractor GWAS: ${SCENARIO} / ${PANEL} panel ==="
  # 1. Convert this panel's FLARE output to Tractor inputs (main-flow steps 11-12).
  python "$(dirname "$0")/../flare/postprocess/11_fill_uncalled_tracts.py" \
      --flare_anc_vcf "${EVAL_LAI_DIR}/${SCENARIO}/${PANEL}.${SIM_CHROM}.anc.vcf.gz" \
      --target_sites  "${EVAL_SIM_DIR}/${SCENARIO}/sim_sites.tsv" \
      --out "anc_${PANEL}.tsv.gz"
  python "$(dirname "$0")/../flare/postprocess/12_flare_to_tractor.py" \
      --phased_vcf "${EVAL_SIM_DIR}/${SCENARIO}/sim_${SCENARIO}.phased.vcf.gz" \
      --ancestry_table "anc_${PANEL}.tsv.gz" \
      --chrom "${SIM_CHROM}" \
      --out_prefix "tractor_in_${PANEL}"

  # 2. Run Tractor's ancestry-aware regression (quantitative trait here).
  python "${TRACTOR_DIR}/RunTractor.py" \
      --hapdose "tractor_in_${PANEL}" \
      --phe ./pheno.tsv \
      --method linear \
      --out "gwas_${PANEL}.tsv"

  gsutil cp "gwas_${PANEL}.tsv" "${EVAL_GWAS_DIR}/${SCENARIO}/"
done
echo "Both panels done. Compare with 05_gwas_power_bias.py"
