#!/usr/bin/env bash
# Step 8 — build the FLARE reference panels from the joint-called
# 1KG-HGDP + 50 MX Biobank WGS callset.
#
# Two panels so you can benchmark panel fit (your accuracy comparison):
#   A) TARGETED   : EUR + AFR + AMR samples that are >95% homogeneous
#                   (sample lists from your LOCAL ADMIXTURE run), MXB included
#                   for a well-matched Amerindigenous AMR anchor.
#   B) COMPARISON : same EUR/AFR, but a broader AMR set WITHOUT the MXB samples
#                   (typical admixed AMR reference) — the "less-fitting" panel.
#
# FLARE wants a phased reference VCF plus a two-column ref sample->panel map.
set -euo pipefail
source "$(dirname "$0")/../../config/config.sh"

WORK="${HOME}/flare_ref"; mkdir -p "${WORK}"; cd "${WORK}"

# 1. Pull the combined callset and the homogeneous sample lists.
gsutil -m cp "${REF_COMBINED_VCF}"* .
COMBINED="$(basename "${REF_COMBINED_VCF}")"
for f in "${REF_HOMOG_EUR_SAMPLES}" "${REF_HOMOG_AFR_SAMPLES}" \
         "${REF_HOMOG_AMR_SAMPLES}" "${REF_COMPARISON_AMR_SAMPLES}"; do
  gsutil cp "$f" .
done
EUR=$(basename "${REF_HOMOG_EUR_SAMPLES}")
AFR=$(basename "${REF_HOMOG_AFR_SAMPLES}")
AMR=$(basename "${REF_HOMOG_AMR_SAMPLES}")
AMR_BROAD=$(basename "${REF_COMPARISON_AMR_SAMPLES}")

# 2. Build the FLARE ref map (sample <TAB> panel-label) for each panel.
build_map () {  # $1=out  $2..=("label" file) pairs
  local out="$1"; shift; : > "$out"
  while (($#)); do local lab="$1" file="$2"; shift 2
    awk -v L="$lab" '{print $1"\t"L}' "$file" >> "$out"; done
}
build_map targeted.map    EUR "$EUR" AFR "$AFR" AMR "$AMR"
build_map comparison.map  EUR "$EUR" AFR "$AFR" AMR "$AMR_BROAD"

# 3. Subset + phase each panel.  (Skip phasing if your combined callset is
#    already phased — check with `bcftools view -h | grep '##phasing'`.)
subset_and_phase () {  # $1=map  $2=out_prefix
  local map="$1" out="$2"
  cut -f1 "$map" | sort -u > "${out}.keep"
  bcftools view -S "${out}.keep" --force-samples "${COMBINED}" -Oz -o "${out}.subset.vcf.gz"
  bcftools index -t "${out}.subset.vcf.gz"
  # Phase with SHAPEIT5 (per-chrom in practice; shown whole-genome for brevity).
  # If already phased, just: mv "${out}.subset.vcf.gz" "${out}.phased.vcf.gz"
  SHAPEIT5_phase_common \
      --input "${out}.subset.vcf.gz" \
      --map   "<${GENETIC_MAP_DIR}/chr@.b38.gmap>" \
      --output "${out}.phased.vcf.gz" \
      --thread "${N_THREADS}" \
      --region "<chrN or leave for a per-chrom loop>"
  bcftools index -t "${out}.phased.vcf.gz"
}
subset_and_phase targeted.map    targeted_panel
subset_and_phase comparison.map  comparison_panel

# 4. Stash both panels.
gsutil -m cp targeted_panel.phased.vcf.gz*   "${FLARE_REF_DIR}/targeted/"
gsutil -m cp targeted.map                    "${FLARE_REF_DIR}/targeted/"
gsutil -m cp comparison_panel.phased.vcf.gz* "${FLARE_REF_DIR}/comparison/"
gsutil -m cp comparison.map                  "${FLARE_REF_DIR}/comparison/"
echo "Reference panels → ${FLARE_REF_DIR}/{targeted,comparison}/"

# NOTE on 2-way vs 3-way: FLARE reads the panel labels present in the ref map.
# For a 2-way AFR-EUR run, pass a map containing only EUR/AFR rows; for 3-way,
# use the full EUR/AFR/AMR map.  Step 10 derives the right map per cohort.
