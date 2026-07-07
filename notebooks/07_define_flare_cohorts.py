"""
Step 7 — finalize the FLARE cohorts by folding qualifying multiracial
participants into the 2-way / 3-way groups using the supervised-ADMIXTURE Q.

Rule of thumb (tune the thresholds to your data / manuscript definitions):
    * 3-way AFR-EUR-AMR  : AFR+EUR+AMR mass >= THREE_WAY_MIN_SUM AND each of the
                           three >= EACH_MIN, and the "other 4" superpops small.
    * 2-way AFR-EUR      : AFR+EUR mass >= TWO_WAY_MIN_SUM AND AMR < AMR_MAX.
    * Participants dominated by SAS/EAS/MEN/OCN are left in "multiracial_other".

Final cohorts written:
    flare_2way_afr_eur.samples.txt
    flare_3way_afr_eur_amr.samples.txt
"""

import os
import pandas as pd

GA_COHORT_DIR = os.environ["GA_COHORT_DIR"]

# ---- tunable thresholds -----------------------------------------------------
THREE_WAY_MIN_SUM = 0.90    # AFR+EUR+AMR must dominate
EACH_MIN          = 0.05    # each of AFR/EUR/AMR at least this
TWO_WAY_MIN_SUM   = 0.90    # AFR+EUR must dominate
AMR_MAX           = 0.05    # ...with negligible AMR to count as 2-way
# ---------------------------------------------------------------------------

os.system(f"gsutil -m cp {GA_COHORT_DIR}/multiracial_admixture_Q.tsv "
          f"{GA_COHORT_DIR}/superpop_column_order.txt "
          f"{GA_COHORT_DIR}/flare_2way.provisional.tsv "
          f"{GA_COHORT_DIR}/flare_3way.provisional.tsv /tmp/")

cols = [l.strip() for l in open("/tmp/superpop_column_order.txt")]  # e.g. AFR,AMR,EAS,EUR,MEN,OCN,SAS
q = pd.read_csv("/tmp/multiracial_admixture_Q.tsv", sep="\t", header=None,
                names=["Sample"] + cols)

afr, eur, amr = q["AFR"], q["EUR"], q["AMR"]
three_way = (afr + eur + amr >= THREE_WAY_MIN_SUM) & (afr >= EACH_MIN) & \
            (eur >= EACH_MIN) & (amr >= EACH_MIN)
two_way   = (afr + eur >= TWO_WAY_MIN_SUM) & (amr < AMR_MAX)

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

pd.Series(final_2way).to_csv("/tmp/flare_2way_afr_eur.samples.txt", index=False, header=False)
pd.Series(final_3way).to_csv("/tmp/flare_3way_afr_eur_amr.samples.txt", index=False, header=False)
print(f"FINAL 2-way AFR-EUR    : {len(final_2way):,}")
print(f"FINAL 3-way AFR-EUR-AMR: {len(final_3way):,}")

os.system(f"gsutil cp /tmp/flare_2way_afr_eur.samples.txt "
          f"/tmp/flare_3way_afr_eur_amr.samples.txt {GA_COHORT_DIR}/")
print(f"Final FLARE cohorts → {GA_COHORT_DIR}/")
