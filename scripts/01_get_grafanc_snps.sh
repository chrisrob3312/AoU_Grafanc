#!/usr/bin/env bash
# Step 1: clone GrafAnc and pull out the ancestry-informative SNP list.
# Run this in an AoU Researcher Workbench "Cloud Analysis" terminal
# (the standard Jupyter env is fine — no Hail/Dataproc needed for this step).
set -euo pipefail

# --- you fill these in -----------------------------------------------------
GRAFANC_REPO="<grafanc_repo_url e.g. https://github.com/jimmy-penn/GrafAnc.git>"
WORKSPACE_NAME="<copy and paste workspace name>"
# ---------------------------------------------------------------------------

# AoU provides $WORKSPACE_BUCKET inside the RW.  Everything we want to keep
# across notebook restarts must live in the bucket, not on the persistent disk.
: "${WORKSPACE_BUCKET:?WORKSPACE_BUCKET not set — are you running inside the RW?}"

WORK=~/grafanc
mkdir -p "$WORK"
cd "$WORK"

# 1. Pull the tool
if [[ ! -d GrafAnc ]]; then
  git clone "$GRAFANC_REPO" GrafAnc
fi
cd GrafAnc

# 2. Locate the SNP definition file.  In NCBI GRAF-pop this was
#    `FpSNPs.txt` shipped with the binary; GrafAnc bundles its own
#    ~280-SNP AIM panel — check the repo and update the path below.
SNP_TABLE="$(find . -maxdepth 3 -iname 'AncSnp*.txt' -o -iname 'FpSNPs*.txt' -o -iname '*aim*snp*.tsv' | head -1)"
if [[ -z "$SNP_TABLE" ]]; then
  echo "Could not auto-locate the AIM SNP file in $GRAFANC_REPO — open the repo and set SNP_TABLE manually." >&2
  exit 1
fi
echo "Using SNP table: $SNP_TABLE"
wc -l "$SNP_TABLE"

# 3. Normalize to a tidy TSV: chrom, pos_b38, rsid, ref, alt
#    (Adjust the awk columns to match whatever shape the file has.)
mkdir -p "$WORK/out"
awk 'BEGIN{OFS="\t"; print "chrom","pos","rsid","ref","alt"}
     NR>1 && $1!~/^#/ {print "chr"$2, $3, $1, $4, $5}' \
     "$SNP_TABLE" > "$WORK/out/grafanc_aim_snps.b38.tsv"

# 4. Stash a copy in the workspace bucket so the Hail step can read it.
gsutil cp "$WORK/out/grafanc_aim_snps.b38.tsv" \
          "$WORKSPACE_BUCKET/grafanc/grafanc_aim_snps.b38.tsv"

echo "Done. SNP list at $WORKSPACE_BUCKET/grafanc/grafanc_aim_snps.b38.tsv"
