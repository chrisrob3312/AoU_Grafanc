#!/usr/bin/env bash
# Step 4 — run GrafAnc on the exported PLINK fileset.
# Run in an RW Cloud Analysis terminal (CPU-only is fine; this is fast).
set -euo pipefail

: "${WORKSPACE_BUCKET:?WORKSPACE_BUCKET not set}"

WORK=~/grafanc
PLINK_PREFIX="$WORK/data/aou_v8_grafanc_aims"
OUT_DIR="$WORK/out"
mkdir -p "$(dirname "$PLINK_PREFIX")" "$OUT_DIR"

# 1. Pull the PLINK fileset down from the workspace bucket.
gsutil -m cp "$WORKSPACE_BUCKET/grafanc/aou_v8_grafanc_aims."{bed,bim,fam} \
        "$(dirname "$PLINK_PREFIX")/"

# 2. Build GrafAnc if not already built.
cd "$WORK/GrafAnc"
if [[ ! -x ./grafanc ]]; then
  make            # or `cmake . && make` — check the repo's README
fi

# 3. Run.  Flags below mirror the GrafAnc README — confirm against the
#    version you cloned, since flag names occasionally change.
./grafanc \
    --bfile  "$PLINK_PREFIX" \
    --out    "$OUT_DIR/aou_v8_grafanc" \
    --threads "<num_threads e.g. 8>"

# 4. Stash results back in the bucket so they outlive the VM.
gsutil -m cp "$OUT_DIR"/aou_v8_grafanc* \
        "$WORKSPACE_BUCKET/grafanc/results/"

echo "Results: $WORKSPACE_BUCKET/grafanc/results/"
