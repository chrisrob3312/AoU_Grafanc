"""
Step 2 — extract the 282,424 GrafAnc ancestry SNPs for EVERY AoU v9 participant
that has imputed-array OR WGS data, and union the two sources.

Run on a Hail / Dataproc cluster in the RW.  Cost is driven by VDS/MT I/O,
not the tiny SNP count.

Strategy
--------
* Pull the ancestry SNPs from the WGS VDS and from the imputed MT separately.
* For participants present in both, prefer WGS (higher quality); imputed fills
  in the array-only participants.
* Emit one dense MatrixTable of ancestry-SNP genotypes over the union cohort.
"""

import os
import hail as hl

WORKSPACE_BUCKET = os.environ["WORKSPACE_BUCKET"]
WGS_VDS_PATH     = os.environ["WGS_VDS_PATH"]
IMPUTED_MT_PATH  = os.environ["IMPUTED_MT_PATH"]
GA_SNP_DIR       = os.environ["GA_SNP_DIR"]
GA_QC_DIR        = os.environ["GA_QC_DIR"]

SNP_TSV = f"{GA_SNP_DIR}/anc_snps.b38.tsv"          # from step 1
MT_WGS  = f"{GA_QC_DIR}/anc_snps_wgs.mt"
MT_IMP  = f"{GA_QC_DIR}/anc_snps_imputed.mt"
MT_UNION= f"{GA_QC_DIR}/anc_snps_union.mt"

hl.init(default_reference="GRCh38", idempotent=True)

# ---- ancestry-SNP intervals -------------------------------------------------
snps = hl.import_table(SNP_TSV, types={"chrom": hl.tstr, "pos": hl.tint32, "rsid": hl.tstr})
snps = snps.annotate(locus=hl.locus(snps.chrom, snps.pos, "GRCh38")).key_by("locus")
intervals = snps.annotate(
    interval=hl.locus_interval(snps.locus.contig, snps.locus.position,
                               snps.locus.position + 1, reference_genome="GRCh38")
).key_by("interval").select()
print(f"Ancestry SNP loci: {snps.count()}")


def extract_from_vds(vds_path):
    vds = hl.vds.read_vds(vds_path)
    vds = hl.vds.filter_intervals(vds, intervals, keep=True)
    mt = hl.vds.to_dense_mt(vds)
    return mt


def extract_from_mt(mt_path):
    mt = hl.read_matrix_table(mt_path)
    return hl.filter_intervals(
        mt,
        [hl.parse_locus_interval(f"{r.chrom}:{r.pos}-{r.pos + 1}", "GRCh38")
         for r in snps.select("chrom", "pos").collect()],
        keep=True,
    )


# ---- WGS --------------------------------------------------------------------
mt_wgs = extract_from_vds(WGS_VDS_PATH)
mt_wgs = mt_wgs.annotate_rows(rsid=snps[mt_wgs.locus].rsid)
mt_wgs = mt_wgs.annotate_cols(src=hl.literal("wgs"))
mt_wgs = mt_wgs.checkpoint(MT_WGS, overwrite=True)
print(f"WGS: {mt_wgs.count_rows()} SNPs × {mt_wgs.count_cols()} samples")

# ---- Imputed array ----------------------------------------------------------
mt_imp = extract_from_mt(IMPUTED_MT_PATH)
mt_imp = mt_imp.annotate_rows(rsid=snps[mt_imp.locus].rsid)
mt_imp = mt_imp.annotate_cols(src=hl.literal("imputed"))
mt_imp = mt_imp.checkpoint(MT_IMP, overwrite=True)
print(f"Imputed: {mt_imp.count_rows()} SNPs × {mt_imp.count_cols()} samples")

# ---- Union of participants, WGS wins where both exist -----------------------
wgs_samples = mt_wgs.cols().s.collect()
mt_imp_only = mt_imp.filter_cols(~hl.literal(set(wgs_samples)).contains(mt_imp.s))
print(f"Imputed-only participants (no WGS): {mt_imp_only.count_cols()}")

# Align entry schemas to just GT before union (drop caller-specific fields).
mt_wgs_g = mt_wgs.select_entries("GT").select_rows("rsid").select_cols("src")
mt_imp_g = mt_imp_only.select_entries("GT").select_rows("rsid").select_cols("src")

mt_union = mt_wgs_g.union_cols(mt_imp_g, row_join_type="outer")
mt_union = mt_union.checkpoint(MT_UNION, overwrite=True)
print(f"UNION cohort: {mt_union.count_rows()} SNPs × {mt_union.count_cols()} participants")
print(f"Wrote {MT_UNION}")
