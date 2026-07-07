#!/usr/bin/env bash
# Step 9 — prepare and phase the AoU target genotypes for a FLARE cohort.
#
# FLARE needs a PHASED target VCF over genome-wide markers (NOT the 282k
# ancestry SNPs — those were only for GrafAnc categorization).  Here we pull
# genome-wide WGS/imputed genotypes for one cohort's participants and phase.
#
# Usage:  ./09_phase_target.sh <cohort_samples.txt> <out_prefix>
#   e.g.  ./09_phase_target.sh flare_3way_afr_eur_amr.samples.txt target_3way
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"

COHORT_KEEP="${1:?pass a cohort sample list, e.g. flare_3way_afr_eur_amr.samples.txt}"
OUT_PREFIX="${2:?pass an output prefix, e.g. target_3way}"

WORK="${HOME}/flare_target/${OUT_PREFIX}"; mkdir -p "${WORK}"; cd "${WORK}"
gsutil cp "${GA_COHORT_DIR}/${COHORT_KEEP}" ./keep.txt

# 1. Extract genome-wide biallelic SNPs for this cohort.  Source: the AoU WGS
#    ACAF-threshold callset (recommended for LAI — common variants, QC'd) or the
#    imputed callset for array-only participants.  Confirm the v9 path.
ACAF_VCF="<gs://.../wgs/acaf_threshold/vcf/chr@.vcf.gz>"   # per-chrom, {@}=chr

for chr in $(seq 1 22); do
  bcftools view -S ./keep.txt --force-samples \
      -v snps -m2 -M2 \
      "${ACAF_VCF/@/chr${chr}}" \
      -Oz -o "chr${chr}.cohort.vcf.gz"
  bcftools index -t "chr${chr}.cohort.vcf.gz"

  # 2. Phase per chromosome (SHAPEIT5).  If AoU already provides phased WGS for
  #    v9, skip this and use those files directly.
  SHAPEIT5_phase_common \
      --input  "chr${chr}.cohort.vcf.gz" \
      --map    "<${GENETIC_MAP_DIR}/chr${chr}.b38.gmap>" \
      --region "chr${chr}" \
      --output "chr${chr}.${OUT_PREFIX}.phased.vcf.gz" \
      --thread "${N_THREADS}"
  bcftools index -t "chr${chr}.${OUT_PREFIX}.phased.vcf.gz"

  gsutil cp "chr${chr}.${OUT_PREFIX}.phased.vcf.gz"* "${FLARE_ROOT}/target/${OUT_PREFIX}/"
done
echo "Phased target → ${FLARE_ROOT}/target/${OUT_PREFIX}/ (per chrom)"
