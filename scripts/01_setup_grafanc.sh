#!/usr/bin/env bash
# Step 1 — clone + build GrafAnc and sanity-check the ancestry-SNP table.
# Run in an RW "Cloud Analysis" terminal (CPU only; fast).
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"

mkdir -p "$(dirname "$GRAFANC_DIR")"

# 1. Clone (or update) GrafAnc.
if [[ ! -d "$GRAFANC_DIR" ]]; then
  git clone "$GRAFANC_REPO" "$GRAFANC_DIR"
else
  git -C "$GRAFANC_DIR" pull --ff-only || true
fi

# 2. Build the C++ executable.
cd "$GRAFANC_DIR/cpp"
make
test -x ./grafanc || { echo "grafanc build failed"; exit 1; }

# GrafAnc needs its data/ dir; set GRAFPATH so it can be called from anywhere.
export GRAFPATH="$GRAFANC_DIR/cpp"
echo "export GRAFPATH=$GRAFPATH" >> "$HOME/.bashrc"

# 3. Confirm the 282,424 ancestry SNPs are present.
#    Columns: chrom, GRCh37 pos, GRCh38 pos, RS ID, ...
if [[ ! -f "$ANC_SNP_TABLE" ]]; then
  echo "AncSnpPopAFs.txt not found at $ANC_SNP_TABLE — check the repo layout." >&2
  exit 1
fi
N_SNP=$(( $(wc -l < "$ANC_SNP_TABLE") - 1 ))
echo "GrafAnc ancestry SNPs available: $N_SNP  (expected ~282424)"

# 4. Emit a tidy GRCh38 SNP list for the Hail extraction step and stash it in
#    the bucket.  chrom | pos_b38 | rsid  (biallelic ancestry SNPs).
mkdir -p "$HOME/tools/out"
awk 'BEGIN{OFS="\t"; print "chrom","pos","rsid"}
     NR>1 && $1!~/^#/ { c=$1; if (c !~ /^chr/) c="chr"c; print c,$3,$4 }' \
     "$ANC_SNP_TABLE" > "$HOME/tools/out/anc_snps.b38.tsv"

gsutil cp "$HOME/tools/out/anc_snps.b38.tsv" "$GA_SNP_DIR/anc_snps.b38.tsv"
echo "SNP list → $GA_SNP_DIR/anc_snps.b38.tsv"
