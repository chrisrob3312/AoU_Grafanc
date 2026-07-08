#!/usr/bin/env bash
# Tier 1 (step 1) — simulate admixed individuals with GROUND-TRUTH local ancestry.
#
# haptools simgenotype draws admixed haplotypes from reference founders under a
# demographic model and writes the TRUE per-haplotype local-ancestry breakpoints
# (.bp).  We then paint these simulated people with each panel (step 2 of the
# main FLARE flow) and score against the truth.
#
# CRITICAL: founders used to simulate must be HELD OUT of the FLARE reference
# panel, or accuracy is inflated.  We split reference samples 50/50 by ID:
# one half = simulation founders, other half = the panel FLARE paints against.
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"
source "$(dirname "$0")/eval_config.sh"

WORK="${HOME}/eval_sim"; mkdir -p "${WORK}"; cd "${WORK}"

# 1. Pull the combined reference callset + homogeneous sample lists.
gsutil -m cp "${REF_COMBINED_VCF}"* .
COMBINED="$(basename "${REF_COMBINED_VCF}")"
for f in "${REF_HOMOG_EUR_SAMPLES}" "${REF_HOMOG_AFR_SAMPLES}" "${REF_HOMOG_AMR_SAMPLES}"; do
  gsutil cp "$f" .
done

# 2. Deterministic 50/50 founder-vs-panel split per ancestry (sort, take alt lines).
#    (No RNG so the split is reproducible across resumes.)
split_founders () {  # $1=samples file  -> ${1}.founders / ${1}.panel
  sort "$1" | awk 'NR%2==1{print > FILENAME".founders"} NR%2==0{print > FILENAME".panel"}' FILENAME="$1"
}
for f in *_homog95.samples.txt; do split_founders "$f"; done
cat ./*_homog95.samples.txt.founders > founders.keep
cat ./*_homog95.samples.txt.panel    > panel.keep

# 3. Subset the founder VCF (what simgenotype samples from).
bcftools view -S founders.keep --force-samples "${COMBINED}" \
  -r "${SIM_CHROM}" -Oz -o founders.${SIM_CHROM}.vcf.gz
bcftools index -t founders.${SIM_CHROM}.vcf.gz

# 4. Build a sample->superpop map for the founders (simgenotype needs population labels).
#    AFR/EUR/AMR from the homogeneous lists.
: > founders.sampleinfo
awk '{print $1"\tAFR"}' afr_homog95.samples.txt.founders >> founders.sampleinfo
awk '{print $1"\tEUR"}' eur_homog95.samples.txt.founders >> founders.sampleinfo
awk '{print $1"\tAMR"}' amr_homog95.samples.txt.founders >> founders.sampleinfo

# 5. Simulate each scenario.  A .dat model file encodes generations + per-pop
#    admixture fractions; we write one per scenario.
for scen in "${SIM_SCENARIOS[@]}"; do
  label="${scen%%:*}"; props="${scen#*:}"
  IFS=',' read -r p_afr p_eur p_amr <<< "${props}"
  cat > model_${label}.dat <<EOF
${SIM_N_INDIV}  Admixed  AFR  EUR  AMR
${SIM_GENERATIONS}  0  ${p_afr}  ${p_eur}  ${p_amr}
EOF

  "${HAPTOOLS_BIN}" simgenotype \
      --model model_${label}.dat \
      --mapdir "<${GENETIC_MAP_DIR}>" \
      --chroms "${SIM_CHROM/chr/}" \
      --ref_vcf founders.${SIM_CHROM}.vcf.gz \
      --sample_info founders.sampleinfo \
      --out sim_${label}.vcf.gz
  # sim_${label}.bp holds the TRUE per-haplotype local ancestry breakpoints.
  bcftools index -t sim_${label}.vcf.gz
  gsutil -m cp sim_${label}.vcf.gz* sim_${label}.bp "${EVAL_SIM_DIR}/${label}/"
done

# 6. Stash the held-out PANEL so FLARE paints against founders it never saw.
bcftools view -S panel.keep --force-samples "${COMBINED}" \
  -Oz -o eval_panel.heldout.vcf.gz
bcftools index -t eval_panel.heldout.vcf.gz
gsutil -m cp eval_panel.heldout.vcf.gz* panel.keep "${EVAL_SIM_DIR}/"
echo "Simulated cohorts + ground-truth .bp → ${EVAL_SIM_DIR}/"
echo "NEXT: phase sim VCFs (step 9) and run FLARE (step 10) with BOTH panels,"
echo "      restricting the ref panel to panel.keep so founders stay held out."
