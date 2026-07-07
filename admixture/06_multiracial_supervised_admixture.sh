#!/usr/bin/env bash
# Step 6 — resolve the GrafAnc "Multiracial" (AncGroupID 800) group.
#
# GrafAnc only models 3 continental components (E/F/A), so AMR admixture is not
# directly readable from its output.  To find which multiracial participants are
# predominantly a 3-way AFR-EUR-AMR mixture, we run SUPERVISED ADMIXTURE with 7
# superpopulation reference sets:  AMR, AFR, EUR, SAS, EAS, MEN, OCN.
#
# Output: per-participant Q proportions across the 7 superpops.  Step 7 then
# pulls participants whose mass sits mostly on {AFR, EUR, AMR} into the 3-way
# FLARE cohort (and, optionally, near-2-way AFR-EUR ones into the 2-way cohort).
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"

WORK="${HOME}/multiracial_admix"
mkdir -p "${WORK}"; cd "${WORK}"

# ---- Inputs you provide -----------------------------------------------------
# Reference-labelled samples for the 7 superpops.  Build this from 1KG+HGDP
# (+ MXB for AMR) with a .pop file where reference samples carry their superpop
# label and the AoU multiracial target samples carry "-" (unknown) so ADMIXTURE
# runs in supervised mode.
REF_MERGED_PLINK="<gs://.../7superpop_ref.{bed,bim,fam} prefix>"   # ref panel
POP_FILE="<gs://.../7superpop.pop>"   # one label per fam row: AMR/AFR/EUR/SAS/EAS/MEN/OCN or "-"
# Multiracial AoU target genotypes at the ancestry SNPs (from step 3 PLINK,
# subset with --keep to multiracial.samples.txt).
MULTIRACIAL_KEEP="${GA_COHORT_DIR}/multiracial.samples.txt"
K=7

# 1. Pull inputs.
gsutil -m cp "${REF_MERGED_PLINK%.*}".{bed,bim,fam} .  2>/dev/null || \
  gsutil -m cp ${REF_MERGED_PLINK}* .
gsutil cp "${POP_FILE}" ./7superpop.pop
gsutil cp "${MULTIRACIAL_KEEP}" ./multiracial.samples.txt
gsutil -m cp "${GA_QC_DIR}/aou_v9_anc_snps_qc."{bed,bim,fam} .

# 2. Merge multiracial target (ancestry SNPs) with the reference panel on the
#    shared ancestry SNPs, keeping only intersecting variants.
plink --bfile aou_v9_anc_snps_qc --keep-fam multiracial.samples.txt \
      --make-bed --out multiracial_targets
plink --bfile 7superpop_ref --bmerge multiracial_targets \
      --geno 0.05 --make-bed --out admix_input
# (Resolve any flip/exclude with the .missnp file if --bmerge complains.)

# 3. Supervised ADMIXTURE (reference labels fixed, targets estimated).
#    Requires admix_input.pop aligned to admix_input.fam row order.
cp 7superpop.pop admix_input.pop
admixture --supervised -j"${N_THREADS}" admix_input.bed ${K}

# admix_input.${K}.Q  = per-sample proportions in ref-panel column order.
# Column order follows the sorted unique labels in the .pop file — record it:
sort -u 7superpop.pop | grep -v '^-$' > superpop_column_order.txt

# 4. Attach sample IDs and stash.
paste <(awk '{print $2}' admix_input.fam) admix_input.${K}.Q > multiracial_admixture_Q.tsv
gsutil cp multiracial_admixture_Q.tsv "${GA_COHORT_DIR}/multiracial_admixture_Q.tsv"
gsutil cp superpop_column_order.txt   "${GA_COHORT_DIR}/superpop_column_order.txt"
echo "Multiracial ADMIXTURE Q → ${GA_COHORT_DIR}/multiracial_admixture_Q.tsv"
