"""
Step 3 — light QC on the ancestry SNPs and export to PLINK for GrafAnc.

Because there is no lab-curated QC-pass list for v9, we do our own quick pass.
GrafAnc explicitly tolerates missingness and needs no LD pruning, so this is
intentionally minimal — quality + call rate + allele sanity only.

We deliberately do NOT apply a global HWE filter: the union cohort is
genetically structured (that's the whole point), so HWE departures reflect
ancestry, not genotyping error.  If you want an HWE screen, run it WITHIN a
homogeneous subgroup, not across the whole cohort.
"""

import os
import hail as hl

GA_QC_DIR = os.environ["GA_QC_DIR"]
MIN_CALL_RATE = float(os.environ.get("QC_MIN_CALL_RATE", 0.95))
WGS_MIN_GQ    = int(os.environ.get("QC_WGS_MIN_GQ", 20))
WGS_MIN_DP    = int(os.environ.get("QC_WGS_MIN_DP", 10))

MT_UNION  = f"{GA_QC_DIR}/anc_snps_union.mt"
PLINK_OUT = f"{GA_QC_DIR}/aou_v9_anc_snps_qc"       # .bed/.bim/.fam prefix

hl.init(default_reference="GRCh38", idempotent=True)
mt = hl.read_matrix_table(MT_UNION)
n0 = mt.count_rows()

# 1. Entry-level quality: only enforce GQ/DP where those fields exist (WGS).
#    Imputed GT carried through as-is (imputation quality was screened upstream
#    in step 2's source, and hard-called array sites are high quality).
if "GQ" in mt.entry and "DP" in mt.entry:
    mt = mt.filter_entries(
        hl.is_missing(mt.GQ) | ((mt.GQ >= WGS_MIN_GQ) & (mt.DP >= WGS_MIN_DP)),
        keep=True,
    )

# 2. Biallelic SNPs only (ancestry SNPs are biallelic by construction; this
#    drops any multiallelic contamination from the join).
mt = mt.filter_rows(hl.len(mt.alleles) == 2)
mt = mt.filter_rows(hl.is_snp(mt.alleles[0], mt.alleles[1]))

# 3. Per-variant call rate across the whole cohort.
mt = hl.variant_qc(mt, name="vqc")
mt = mt.filter_rows(mt.vqc.call_rate >= MIN_CALL_RATE)

# 4. Drop monomorphic sites (uninformative, and can confuse PLINK export).
mt = mt.filter_rows((mt.vqc.AF[1] > 0) & (mt.vqc.AF[1] < 1))

n1 = mt.count_rows()
print(f"Ancestry SNPs: {n0} → {n1} after QC "
      f"(call_rate≥{MIN_CALL_RATE}, biallelic SNP, polymorphic)")

# 5. Per-sample call rate report (informational — GrafAnc handles missingness,
#    but flag samples with very few genotyped ancestry SNPs; <100 → no result).
mt = hl.sample_qc(mt, name="sqc")
lowcov = mt.aggregate_cols(hl.agg.count_where(mt.sqc.n_called < 10000))
print(f"Participants with <10,000 genotyped ancestry SNPs (lower reliability): {lowcov}")

# 6. Export PLINK for GrafAnc.  varid = rsid so GrafAnc matches by RS ID.
hl.export_plink(mt, PLINK_OUT, ind_id=mt.s, varid=mt.rsid)
print(f"Wrote {PLINK_OUT}.{{bed,bim,fam}}")
