"""
Step 2 — subset the AoU v8 WGS VDS to the GrafAnc AIM SNPs.

Run this in a Hail / Dataproc notebook on the Researcher Workbench.
Recommended cluster: 1 main + <num_workers> preemptible workers, n1-highmem-8.

The VDS has ~245k participants in v8, so a list of ~280 SNPs is tiny;
the cost driver is the I/O over the VDS, not the SNP count.
"""

import os
import hail as hl

# ---- workspace-provided env vars (don't hard-code paths) ------------------
WORKSPACE_BUCKET = os.environ["WORKSPACE_BUCKET"]
WGS_VDS_PATH     = os.environ["WGS_VDS_PATH"]          # AoU sets this for you
# ---------------------------------------------------------------------------

OUT_PREFIX = f"{WORKSPACE_BUCKET}/grafanc"
SNP_TSV    = f"{OUT_PREFIX}/grafanc_aim_snps.b38.tsv"  # written in step 1
MT_OUT     = f"{OUT_PREFIX}/aou_v8_grafanc_aims.mt"

hl.init(default_reference="GRCh38", idempotent=True)

# 1. Read the AIM SNP list and turn it into a (locus, alleles) keyed table.
snps = hl.import_table(
    SNP_TSV,
    types={"chrom": hl.tstr, "pos": hl.tint32,
           "rsid": hl.tstr,  "ref": hl.tstr, "alt": hl.tstr},
)
snps = snps.annotate(
    locus   = hl.locus(snps.chrom, snps.pos, reference_genome="GRCh38"),
    alleles = hl.array([snps.ref, snps.alt]),
)
snps = snps.key_by("locus", "alleles")
print(f"AIM SNPs in panel: {snps.count()}")

# 2. Load the VDS and densify only the rows we care about.  filter_intervals
#    on the variant data is the cheap way in; we then refine by alleles.
vds = hl.vds.read_vds(WGS_VDS_PATH)

intervals = snps.select().distinct()
intervals = intervals.annotate(
    interval = hl.locus_interval(intervals.locus.contig,
                                 intervals.locus.position,
                                 intervals.locus.position + 1,
                                 reference_genome="GRCh38"),
).key_by("interval")

vds = hl.vds.filter_intervals(vds, intervals, keep=True)
mt  = hl.vds.to_dense_mt(vds)

# 3. Keep only the exact (locus, alleles) pairs in the panel.
mt = mt.semi_join_rows(snps)
mt = mt.annotate_rows(rsid = snps[mt.row_key].rsid)

# 4. Write a checkpoint MT we can reuse cheaply downstream.
mt = mt.checkpoint(MT_OUT, overwrite=True)
print(f"Variants kept: {mt.count_rows()}   Samples: {mt.count_cols()}")
print(f"Wrote {MT_OUT}")
