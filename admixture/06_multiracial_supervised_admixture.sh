#!/usr/bin/env bash
# Step 6 — resolve the GrafAnc "Multiracial" (AncGroupID 800) group.
#
# GrafAnc only models 3 continental components (E/F/A), so AMR admixture is not
# directly readable from its output.  To find which multiracial participants are
# predominantly a 3-way AFR-EUR-AMR mixture, we run SUPERVISED ADMIXTURE with 7
# superpopulation reference sets:  AMR, AFR, EUR, SAS, EAS, MEN, OCN.
#
# Marker set (flag, $1):
#   ancsnp       Use the ~282k GrafAnc ancestry SNPs (step 3 QC'd PLINK).
#                Fast, LD already accounted for → no pruning needed.  Fine for
#                global proportions.
#   genomewide   Use a denser, LD-pruned genome-wide common-variant set.  More
#                markers → tighter separation of the "other" superpops
#                (SAS/EAS/MEN/OCN), which sharpens the OTHER_MAX gate in step 7.
#                Slower; needs a genome-wide reference panel + LD pruning.
#   both         Run both; the genome-wide Q becomes the canonical input to
#                step 7 (override with PRIMARY_MARKERSET below).
#
# Output: per-participant Q proportions across the 7 superpops, written per
# marker set AND copied to the canonical file step 7 reads.
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"

MARKERSET="${1:-ancsnp}"            # ancsnp | genomewide | both
PRIMARY_MARKERSET="${PRIMARY_MARKERSET:-genomewide}"   # which Q step 7 consumes when both
K=7

WORK="${HOME}/multiracial_admix"; mkdir -p "${WORK}"; cd "${WORK}"

# ---- Inputs you provide -----------------------------------------------------
# Reference-labelled samples for the 7 superpops.  Supply a KEYED label file
# (two columns: sample_id <TAB> superpop) so we can rebuild the .pop in the exact
# row order of whatever merged fam PLINK produces — positional .pop files break
# if PLINK reorders samples on merge.  Reference rows carry AMR/AFR/EUR/SAS/EAS/
# MEN/OCN; AoU multiracial targets are not listed and become "-" (unknown), which
# is what makes ADMIXTURE supervised (and fixes cluster identity).
#
# ancsnp path:     reference genotypes restricted to the ancestry SNPs.
# genomewide path: reference genotypes genome-wide (from your combined
#                  1KG-HGDP+MXB callset, REF_COMBINED_VCF), pruned to match.
REF_ANCSNP_PLINK="<gs://.../7superpop_ref.ancsnp.{bed,bim,fam} prefix>"
REF_GW_PLINK="<gs://.../7superpop_ref.genomewide.{bed,bim,fam} prefix>"
REF_LABELS="<gs://.../7superpop_ref_labels.tsv>"   # sample_id <TAB> superpop
MULTIRACIAL_KEEP="${GA_COHORT_DIR}/multiracial.samples.txt"

# Genome-wide AoU target source (common, QC'd variants) — the ACAF-threshold
# WGS callset (per-chrom; {@}=chr) and/or the imputed callset.  Confirm v9 path.
ACAF_VCF="<gs://.../wgs/acaf_threshold/vcf/chr@.vcf.gz>"

# LD-pruning params for the genome-wide path (ADMIXTURE wants ~LE markers).
LD_WINDOW=50; LD_STEP=10; LD_R2=0.1; GW_MAF=0.01

gsutil cp "${REF_LABELS}" ./ref_labels.tsv
gsutil cp "${MULTIRACIAL_KEEP}" ./multiracial.samples.txt

# Column order = sorted unique reference labels (ADMIXTURE Q column convention).
cut -f2 ref_labels.tsv | sort -u > superpop_column_order.txt

# =============================================================================
run_admixture () {   # $1 = markerset label
  local ms="$1"
  echo "=== supervised ADMIXTURE (${ms}) ==="

  if [[ "${ms}" == "ancsnp" ]]; then
    # Target = multiracial subset of the QC'd ancestry-SNP PLINK.
    gsutil -m cp "${GA_QC_DIR}/aou_v9_anc_snps_qc."{bed,bim,fam} .
    gsutil -m cp "${REF_ANCSNP_PLINK%.*}".{bed,bim,fam} . 2>/dev/null || \
      gsutil -m cp ${REF_ANCSNP_PLINK}* .
    plink --bfile aou_v9_anc_snps_qc --keep-fam multiracial.samples.txt \
          --make-bed --out targets_${ms}
    local REF_PREFIX
    REF_PREFIX="$(basename "${REF_ANCSNP_PLINK%.*}")"
    # No LD pruning — ancestry SNPs already selected for near-LE.
    plink --bfile "${REF_PREFIX}" --bmerge targets_${ms} \
          --geno 0.05 --make-bed --out admix_${ms}

  elif [[ "${ms}" == "genomewide" ]]; then
    # Target = genome-wide common biallelic SNPs for the multiracial cohort.
    for chr in $(seq 1 22); do
      bcftools view -S multiracial.samples.txt --force-samples \
          -v snps -m2 -M2 -q ${GW_MAF}:minor \
          "${ACAF_VCF/@/chr${chr}}" -Ob -o tgt_chr${chr}.bcf
      bcftools index tgt_chr${chr}.bcf
    done
    bcftools concat -Ob -o targets_gw.bcf tgt_chr*.bcf
    plink --bcf targets_gw.bcf --make-bed --out targets_${ms}
    # Reference genome-wide, then merge on shared sites.
    gsutil -m cp "${REF_GW_PLINK%.*}".{bed,bim,fam} . 2>/dev/null || \
      gsutil -m cp ${REF_GW_PLINK}* .
    local REF_PREFIX
    REF_PREFIX="$(basename "${REF_GW_PLINK%.*}")"
    plink --bfile "${REF_PREFIX}" --bmerge targets_${ms} \
          --geno 0.05 --maf ${GW_MAF} --make-bed --out admix_merged_${ms}
    # LD-prune the MERGED set so ADMIXTURE markers are ~independent.
    plink --bfile admix_merged_${ms} \
          --indep-pairwise ${LD_WINDOW} ${LD_STEP} ${LD_R2} --out prune_${ms}
    plink --bfile admix_merged_${ms} --extract prune_${ms}.prune.in \
          --make-bed --out admix_${ms}
  else
    echo "unknown markerset '${ms}'"; return 1
  fi

  # Supervised ADMIXTURE — the .pop file must align line-for-line with the merged
  # fam row order.  Rebuild it by looking up each fam IID (col 2) in ref_labels:
  # a reference sample gets its superpop label, an AoU target gets "-" (unknown).
  awk 'NR==FNR{lab[$1]=$2; next} {print ($2 in lab) ? lab[$2] : "-"}' \
      ref_labels.tsv admix_${ms}.fam > admix_${ms}.pop

  admixture --supervised -j"${N_THREADS}" admix_${ms}.bed ${K}

  paste <(awk '{print $2}' admix_${ms}.fam) admix_${ms}.${K}.Q \
        > multiracial_admixture_Q.${ms}.tsv
  gsutil cp multiracial_admixture_Q.${ms}.tsv \
        "${GA_COHORT_DIR}/multiracial_admixture_Q.${ms}.tsv"
  echo "  → ${GA_COHORT_DIR}/multiracial_admixture_Q.${ms}.tsv"
}
# =============================================================================

case "${MARKERSET}" in
  ancsnp)     run_admixture ancsnp;     CANON=ancsnp ;;
  genomewide) run_admixture genomewide; CANON=genomewide ;;
  both)       run_admixture ancsnp; run_admixture genomewide; CANON="${PRIMARY_MARKERSET}" ;;
  *) echo "usage: $0 <ancsnp|genomewide|both>"; exit 1 ;;
esac

# Canonical file consumed by step 7 (header-matched by column order).
gsutil cp "${GA_COHORT_DIR}/multiracial_admixture_Q.${CANON}.tsv" \
          "${GA_COHORT_DIR}/multiracial_admixture_Q.tsv"
gsutil cp superpop_column_order.txt "${GA_COHORT_DIR}/superpop_column_order.txt"
echo "Canonical Q for step 7 = ${CANON} → ${GA_COHORT_DIR}/multiracial_admixture_Q.tsv"
