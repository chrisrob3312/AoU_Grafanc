#!/usr/bin/env bash
# Build FLARE/RFMix reference-panel keep-lists and label files from the
# homogeneous per-super-pop sample lists in reference/homog_by_super_pop/.
#
# The homogeneity panel has 5 super-pops (AFR 634, AMR 88, EAS 667, EUR 620,
# SAS 48). The AMR list of 88 INCLUDES the 50 MX Biobank (MXB_*) samples plus 38
# HGDP/1KG Amerindigenous. So the MXB contribution is isolated exactly by
# dropping the MXB_ IDs -> that is the amr_small_homog arm.
#
# (These are de-identified sample IDs / metadata only -- no sequences. They let
# anyone with MX Biobank DUA access reproduce the exact panel composition.)
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
SRC="${HERE}/homog_by_super_pop"
OUT="${HERE}/panels"; mkdir -p "${OUT}"

# 1. Split AMR into MXB vs non-MXB.
grep    '^MXB_' "${SRC}/AMR_ids.txt" > "${OUT}/AMR_mxb.txt"
grep -v '^MXB_' "${SRC}/AMR_ids.txt" > "${OUT}/AMR_nonmxb.txt"

# 2. Keyed sample->superpop label file (RFMix classes / FLARE ref map /
#    supervised ADMIXTURE all consume this). 5 super-pops.
: > "${OUT}/ref_labels.tsv"
for pop in AFR AMR EAS EUR SAS; do
  awk -v P="$pop" '{print $1"\t"P}' "${SRC}/${pop}_ids.txt" >> "${OUT}/ref_labels.tsv"
done

# 3. Panel keep-lists for the LAI arms. The AoU Latin American cohort is painted
#    with AFR/EUR/AMR, so each arm = AFR + EUR + a choice of AMR reference.
#    (EAS/SAS stay in ref_labels for other cohorts / completeness, but are not in
#    the AFR-EUR-AMR painting arms.)
build_arm () {  # $1=arm name  $2=AMR source file
  cat "${SRC}/AFR_ids.txt" "${SRC}/EUR_ids.txt" "$2" | sort -u > "${OUT}/$1.keep"
  echo "  $1.keep  ($(wc -l < "${OUT}/$1.keep") samples)"
}
echo "Panel arms:"
build_arm mxb_expanded    "${SRC}/AMR_ids.txt"      # AFR+EUR+AMR(all 88, incl MXB)  <- proposal
build_arm amr_small_homog "${OUT}/AMR_nonmxb.txt"   # AFR+EUR+AMR(38, MXB removed)   <- matching, low N

# The other two arms need EXTERNAL AMR references not in the homogeneity panel:
#   amr_large_admixed : full admixed 1KG AMR (MXL/PEL/CLM/PUR) -- size, contaminated
#   aou_default_1kg   : 1KG-style baseline ~ what AoU CDR v9 will ship
# Provide those keep-lists and uncomment:
# cat "${SRC}/AFR_ids.txt" "${SRC}/EUR_ids.txt" "<admixed_amr.txt>" | sort -u > "${OUT}/amr_large_admixed.keep"
# cp "<aou_default_1kg.keep>" "${OUT}/aou_default_1kg.keep"

echo "Wrote ref_labels.tsv + AMR_mxb/AMR_nonmxb splits + arm keep-lists to ${OUT}/"
echo "Counts: AMR total $(wc -l < ${SRC}/AMR_ids.txt), MXB $(wc -l < ${OUT}/AMR_mxb.txt), non-MXB $(wc -l < ${OUT}/AMR_nonmxb.txt)"
