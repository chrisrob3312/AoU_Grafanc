"""
Step 5 — read the GrafAnc output and split participants into the cohorts we
will carry into FLARE.

Cohorts (by AncGroupID):
    107  African American          -> 2-way AFR-EUR candidate
    601  Latin American 1 (LA1)     -> 3-way AFR-EUR-AMR candidate
    602  Latin American 2 (LA2)     -> 3-way AFR-EUR-AMR candidate
    603  Native American            -> AMR-anchored (mostly for the ref/eval side)
    800  Multiracial                -> UNRESOLVED; sent to the supervised
                                       7-superpop ADMIXTURE screen (step 6) to
                                       find those that are truly AFR-EUR-AMR
                                       3-way and pull them into the 3-way cohort.

Everything here is bookkeeping on a single results table, so it runs fine in a
plain Jupyter (pandas) env — no Hail needed.
"""

import os
import pandas as pd

GA_RESULTS_DIR = os.environ["GA_RESULTS_DIR"]
GA_COHORT_DIR  = os.environ["GA_COHORT_DIR"]

RESULTS = f"{GA_RESULTS_DIR}/aou_v9_grafanc_pops.txt"

df = pd.read_csv(RESULTS, sep="\t")
df.columns = [c.strip().lstrip("#") for c in df.columns]   # tidy "#SNPs" etc.
print(f"GrafAnc rows: {len(df):,}")
print(df["AncGroupID"].value_counts().sort_index())

# Reliability flag: GrafAnc gives no continental call under 100 SNPs and is
# less reliable below ~10k.  Keep the flag with the cohort for later filtering.
df["reliable_continental"] = df["SNPs"] >= 10000

COHORTS = {
    "african_american": [107],
    "latin_american_1": [601],
    "latin_american_2": [602],
    "native_american":  [603],
    "multiracial":      [800],
}

os.makedirs("/tmp/cohorts", exist_ok=True)
for name, ids in COHORTS.items():
    sub = df[df["AncGroupID"].isin(ids)]
    # Sample-ID list (for PLINK/VCF --keep) and a full annotated table.
    ids_path = f"/tmp/cohorts/{name}.samples.txt"
    tab_path = f"/tmp/cohorts/{name}.grafanc.tsv"
    sub[["Sample"]].to_csv(ids_path, sep="\t", index=False, header=False)
    sub.to_csv(tab_path, sep="\t", index=False)
    print(f"  {name:20s} n={len(sub):>8,}  reliable={sub['reliable_continental'].sum():>8,}")

# Provisional FLARE assignments (multiracial refined in step 6/7).
#   2-way AFR-EUR : African American
#   3-way AFR-EUR-AMR : LA1 + LA2  (+ qualifying multiracial, added in step 7)
pd.concat([df[df.AncGroupID == 107].assign(flare_mode="2way_afr_eur")]) \
  .to_csv("/tmp/cohorts/flare_2way.provisional.tsv", sep="\t", index=False)
df[df.AncGroupID.isin([601, 602])].assign(flare_mode="3way_afr_eur_amr") \
  .to_csv("/tmp/cohorts/flare_3way.provisional.tsv", sep="\t", index=False)

os.system(f"gsutil -m cp /tmp/cohorts/* {GA_COHORT_DIR}/")
print(f"Cohort lists → {GA_COHORT_DIR}/")
