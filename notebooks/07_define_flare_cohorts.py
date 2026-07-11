"""
Step 7 — finalize the FLARE cohorts by folding qualifying multiracial
participants into the 2-way / 3-way groups using the supervised-ADMIXTURE Q.

Because step 6 runs ADMIXTURE in SUPERVISED mode (reference samples carry fixed
superpop labels), each Q column already corresponds to a known superpopulation —
we don't have to guess which cluster is which. The column order is read from
`superpop_column_order.txt`, so AFR/EUR/AMR/etc. are addressed by name below.

Rule of thumb (tune the thresholds to your data / manuscript definitions):
    * 3-way AFR-EUR-AMR : AFR+EUR+AMR >= THREE_WAY_MIN_SUM, each of the three
                          >= EACH_MIN, AND the OTHER superpops (SAS/EAS/MEN/OCN)
                          stay small — total <= OTHER_MAX and none > OTHER_EACH_MAX.
    * 2-way AFR-EUR     : AFR+EUR >= TWO_WAY_MIN_SUM, AMR < AMR_MAX, and the same
                          OTHER-superpop caps.
    * Everyone else is left in "multiracial_other" (dominated by a non-AFR/EUR/AMR
      superpop, so not a clean AFR-EUR(-AMR) admixture).

Final cohorts written:
    flare_2way_afr_eur.samples.txt
    flare_3way_afr_eur_amr.samples.txt
    flare_admixed_all.samples.txt   (union — for a single combined FLARE run)
"""

import os
import pandas as pd

GA_COHORT_DIR = os.environ["GA_COHORT_DIR"]

# ---- tunable thresholds -----------------------------------------------------
# Target = the ancestries we model in FLARE (AFR, EUR, AMR).
# "Other" = every other superpop in the 7-way ADMIXTURE (SAS, EAS, MEN, OCN):
# these are the components that must be BELOW a threshold for a participant to
# count as a clean AFR-EUR(-AMR) mixture worth painting.
THREE_WAY_MIN_SUM = 0.90    # AFR+EUR+AMR must dominate
EACH_MIN          = 0.05    # each of AFR/EUR/AMR at least this (true 3-way)
TWO_WAY_MIN_SUM   = 0.90    # AFR+EUR must dominate
AMR_MAX           = 0.05    # ...with negligible AMR to count as 2-way
OTHER_MAX         = 0.10    # summed SAS+EAS+MEN+OCN must be <= this
OTHER_EACH_MAX    = 0.05    # and no single "other" superpop above this
# ---------------------------------------------------------------------------

TARGET_ANCS = ["AFR", "EUR", "AMR"]
# ---------------------------------------------------------------------------

os.system(f"gsutil -m cp {GA_COHORT_DIR}/multiracial_admixture_Q.tsv "
          f"{GA_COHORT_DIR}/superpop_column_order.txt "
          f"{GA_COHORT_DIR}/flare_2way.provisional.tsv "
          f"{GA_COHORT_DIR}/flare_3way.provisional.tsv /tmp/")

cols = [l.strip() for l in open("/tmp/superpop_column_order.txt")]  # e.g. AFR,AMR,EAS,EUR,MEN,OCN,SAS
q = pd.read_csv("/tmp/multiracial_admixture_Q.tsv", sep="\t", header=None,
                names=["Sample"] + cols)

afr, eur, amr = q["AFR"], q["EUR"], q["AMR"]

# "Other" superpops = whatever ADMIXTURE columns are not AFR/EUR/AMR.
other_cols = [c for c in cols if c not in TARGET_ANCS]
other_sum  = q[other_cols].sum(axis=1)
other_max_each = q[other_cols].max(axis=1)
other_ok = (other_sum <= OTHER_MAX) & (other_max_each <= OTHER_EACH_MAX)
print(f"'Other' superpops screened: {other_cols}")

three_way = other_ok & (afr + eur + amr >= THREE_WAY_MIN_SUM) & \
            (afr >= EACH_MIN) & (eur >= EACH_MIN) & (amr >= EACH_MIN)
two_way   = other_ok & (afr + eur >= TWO_WAY_MIN_SUM) & (amr < AMR_MAX)

mr_3way = q.loc[three_way, "Sample"]
mr_2way = q.loc[two_way & ~three_way, "Sample"]
print(f"Multiracial → 3-way AFR-EUR-AMR: {len(mr_3way):,}")
print(f"Multiracial → 2-way AFR-EUR   : {len(mr_2way):,}")
print(f"Multiracial → other/unused    : {len(q) - len(mr_3way) - len(mr_2way):,}")

# Base cohorts from GrafAnc (step 5): AA -> 2way; LA1+LA2 -> 3way.
base_2way = pd.read_csv("/tmp/flare_2way.provisional.tsv", sep="\t")["Sample"]
base_3way = pd.read_csv("/tmp/flare_3way.provisional.tsv", sep="\t")["Sample"]

final_2way = pd.Index(base_2way).union(pd.Index(mr_2way))
final_3way = pd.Index(base_3way).union(pd.Index(mr_3way))
# A participant qualifying for both defaults to 3-way (more general model).
final_2way = final_2way.difference(final_3way)

# Combined cohort: everyone we intend to paint.  FLARE can run the 2-way and
# 3-way individuals TOGETHER against the 3-way (EUR/AFR/AMR) reference panel —
# a genuinely 2-way AFR-EUR person simply gets ~0 AMR assigned, which is
# correct.  This is the recommended default (one run, one output, no cohort
# split); the separate 2-way panel run is kept only for the panel comparison
# (step 13) and for anyone who wants a strict no-AMR model.
final_all = pd.Index(final_2way).union(pd.Index(final_3way))

pd.Series(final_2way).to_csv("/tmp/flare_2way_afr_eur.samples.txt", index=False, header=False)
pd.Series(final_3way).to_csv("/tmp/flare_3way_afr_eur_amr.samples.txt", index=False, header=False)
pd.Series(final_all).to_csv("/tmp/flare_admixed_all.samples.txt", index=False, header=False)
print(f"FINAL 2-way AFR-EUR    : {len(final_2way):,}")
print(f"FINAL 3-way AFR-EUR-AMR: {len(final_3way):,}")
print(f"FINAL combined (all)   : {len(final_all):,}  ← recommended single FLARE run vs 3-way panel")

os.system(f"gsutil cp /tmp/flare_2way_afr_eur.samples.txt "
          f"/tmp/flare_3way_afr_eur_amr.samples.txt "
          f"/tmp/flare_admixed_all.samples.txt {GA_COHORT_DIR}/")
print(f"Final FLARE cohorts → {GA_COHORT_DIR}/")
