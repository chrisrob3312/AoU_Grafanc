#!/usr/bin/env bash
# Step 10 — run FLARE local-ancestry inference, per chromosome, for a cohort
# against a chosen reference panel.
#
# Runs the 2-way (AFR-EUR) and 3-way (AFR-EUR-AMR) analyses, each against BOTH
# the targeted and comparison panels so you can benchmark panel fit (step 13).
#
# FLARE: https://github.com/browning-lab/flare
#   java -jar flare.jar ref=<ref.vcf.gz> ref-panel=<map> gt=<target.vcf.gz>
#        map=<plink.map> out=<prefix> nthreads=<n>
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"

FLARE_JAR="${HOME}/tools/flare.jar"
[[ -f "${FLARE_JAR}" ]] || { mkdir -p "$(dirname "$FLARE_JAR")"; wget -O "${FLARE_JAR}" "${FLARE_JAR_URL}"; }

# ---- pick the run: mode × panel (targeted|comparison) -----------------------
# MODE:
#   combined  RECOMMENDED. Paint all admixed individuals (2-way + 3-way together)
#             against the 3-way EUR/AFR/AMR panel. A truly 2-way AFR-EUR person
#             just gets ~0 AMR — correct — so no cohort split is needed.
#   3way      Only the 3-way cohort, against the EUR/AFR/AMR panel.
#   2way      Only the 2-way cohort, against an EUR/AFR-only panel (strict
#             no-AMR model; mainly for the panel comparison in step 13).
MODE="${1:?usage: 10_run_flare.sh <combined|3way|2way> <targeted|comparison>}"
PANEL="${2:?usage: 10_run_flare.sh <combined|3way|2way> <targeted|comparison>}"

case "${MODE}" in
  combined) TGT_PREFIX="target_all";  PANEL_LABELS="EUR,AFR,AMR" ;;
  3way)     TGT_PREFIX="target_3way"; PANEL_LABELS="EUR,AFR,AMR" ;;
  2way)     TGT_PREFIX="target_2way"; PANEL_LABELS="EUR,AFR" ;;
  *) echo "MODE must be combined, 3way, or 2way"; exit 1 ;;
esac

WORK="${HOME}/flare_run/${MODE}_${PANEL}"; mkdir -p "${WORK}"; cd "${WORK}"

# 1. Pull the reference panel and build a mode-specific ref map (subset labels).
gsutil -m cp "${FLARE_REF_DIR}/${PANEL}/"* .
REF_VCF="$(ls *phased.vcf.gz | head -1)"
awk -F'\t' -v keep="${PANEL_LABELS}" 'BEGIN{split(keep,a,","); for(i in a)K[a[i]]=1}
     K[$2]{print}' ./*.map > ref_panel_${MODE}.map

# 2. Run FLARE per chromosome.
mkdir -p out
for chr in $(seq 1 22); do
  gsutil cp "${FLARE_ROOT}/target/${TGT_PREFIX}/chr${chr}.${TGT_PREFIX}.phased.vcf.gz"* . || {
    echo "missing phased target for chr${chr}; run step 9 first"; exit 1; }

  java -Xmx"<flare_heap e.g. 48g>" -jar "${FLARE_JAR}" \
      ref="${REF_VCF}" \
      ref-panel=ref_panel_${MODE}.map \
      gt="chr${chr}.${TGT_PREFIX}.phased.vcf.gz" \
      map="<${GENETIC_MAP_DIR}/chr${chr}.plink.map>" \
      out="out/chr${chr}.${MODE}.${PANEL}" \
      probs=true \
      nthreads="${N_THREADS}"
done

# 3. Stash FLARE output (.anc.vcf.gz + .global.anc.gz + .model + .log).
gsutil -m cp out/* "${FLARE_OUT_DIR}/${MODE}_${PANEL}/"
echo "FLARE ${MODE}/${PANEL} → ${FLARE_OUT_DIR}/${MODE}_${PANEL}/"
