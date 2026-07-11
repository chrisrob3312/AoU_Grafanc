#!/usr/bin/env bash
# Step 4 — run GrafAnc on the QC'd ancestry-SNP PLINK set.
# Run in an RW "Cloud Analysis" terminal.  Fast and CPU-only.
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"
export GRAFPATH="${GRAFANC_DIR}/cpp"

WORK="${HOME}/grafanc_run"
mkdir -p "${WORK}/data" "${WORK}/out"

# 1. Pull the QC'd PLINK set from the bucket.
gsutil -m cp "${GA_QC_DIR}/aou_v9_anc_snps_qc."{bed,bim,fam} "${WORK}/data/"

# 2. Run GrafAnc.  With ~245k+ WGS + array participants this is a large N;
#    raise --maxmem and let GrafAnc auto-batch, or pin --samples per batch.
"${GRAFANC_DIR}/cpp/grafanc" \
    "${WORK}/data/aou_v9_anc_snps_qc.bed" \
    "${WORK}/out/aou_v9_grafanc_pops.txt" \
    --maxmem  "<max_mem_MB e.g. 32000>" \
    --threads "${N_THREADS}"

# 3. Stash results in the bucket.
gsutil -m cp "${WORK}/out/aou_v9_grafanc_pops.txt" "${GA_RESULTS_DIR}/"
echo "GrafAnc results → ${GA_RESULTS_DIR}/aou_v9_grafanc_pops.txt"

# The output table has one row per participant with columns:
#   Sample #SNPs GD1 GD2 GD3 EA1..EA4 AF1..AF3 EU1..EU3 SA1 SA2
#   IC1..IC3 Pe Pf Pa RawPe RawPf RawPa AncGroupID
