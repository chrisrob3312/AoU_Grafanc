"""
Step 3 — intersect the AIM-SNP MT with the lab's v8 QC-pass variant list
and export to PLINK (.bed/.bim/.fam) for GrafAnc.

GrafAnc consumes PLINK binary filesets, so this is the format we want.
"""

import os
import hail as hl

WORKSPACE_BUCKET = os.environ["WORKSPACE_BUCKET"]

# ---- you fill these in ----------------------------------------------------
# A TSV/CSV your lab member maintains.  Expected columns (rename as needed):
#   chrom (e.g. "chr1"), pos (int, GRCh38), ref, alt
# If the list is keyed by rsid only, see the rsid branch below.
QC_PASS_PATH = "<path/to/labmate/qc_pass_variants.tsv>"   # gs:// or local
QC_KEY_TYPE  = "<'locus_alleles' or 'rsid'>"              # which key the list uses
# ---------------------------------------------------------------------------

MT_IN      = f"{WORKSPACE_BUCKET}/grafanc/aou_v8_grafanc_aims.mt"
PLINK_OUT  = f"{WORKSPACE_BUCKET}/grafanc/aou_v8_grafanc_aims"  # no extension

hl.init(default_reference="GRCh38", idempotent=True)
mt = hl.read_matrix_table(MT_IN)
print(f"Start: {mt.count_rows()} variants × {mt.count_cols()} samples")

qc = hl.import_table(QC_PASS_PATH, force_bgz=QC_PASS_PATH.endswith(".bgz"))

if QC_KEY_TYPE == "locus_alleles":
    qc = qc.annotate(
        locus   = hl.locus(qc.chrom, hl.int32(qc.pos), reference_genome="GRCh38"),
        alleles = hl.array([qc.ref, qc.alt]),
    ).key_by("locus", "alleles")
    mt = mt.semi_join_rows(qc)
elif QC_KEY_TYPE == "rsid":
    qc = qc.key_by("rsid")
    mt = mt.filter_rows(hl.is_defined(qc[mt.rsid]))
else:
    raise ValueError(f"Set QC_KEY_TYPE to 'locus_alleles' or 'rsid' (got {QC_KEY_TYPE!r})")

print(f"After QC-pass intersect: {mt.count_rows()} variants")

# Drop monomorphic / all-missing rows that survive the join (rare but possible).
mt = hl.variant_qc(mt, name="vqc")
mt = mt.filter_rows((mt.vqc.AF[1] > 0) & (mt.vqc.AF[1] < 1) & (mt.vqc.call_rate > 0.95))
print(f"After variant_qc filter: {mt.count_rows()} variants")

# PLINK export.  GrafAnc reads sample IDs from the FAM, so use the AoU
# `person_id`-style column key directly.
hl.export_plink(
    mt,
    PLINK_OUT,
    ind_id   = mt.s,           # AoU sample id
    varid    = mt.rsid,        # so GrafAnc can match by rsid
)
print(f"Wrote {PLINK_OUT}.{{bed,bim,fam}}")
