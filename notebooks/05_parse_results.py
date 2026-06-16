"""
Step 5 — load the GrafAnc ancestry calls and join them to the AoU
person table for downstream analysis.
"""

import os
import pandas as pd

WORKSPACE_BUCKET = os.environ["WORKSPACE_BUCKET"]
WORKSPACE_CDR    = os.environ["WORKSPACE_CDR"]  # BigQuery dataset

RESULTS = f"{WORKSPACE_BUCKET}/grafanc/results/aou_v8_grafanc.<grafanc_output_ext e.g. .anc.txt>"

anc = pd.read_csv(RESULTS, sep="\t")
print(anc.head())
print(anc["<ancestry column name e.g. AncestryGroup>"].value_counts())

# Optional: pull demographics from the CDR and merge.
from google.cloud import bigquery
bq = bigquery.Client()
demo = bq.query(f"""
    SELECT person_id,
           gender_concept_id,
           race_concept_id,
           ethnicity_concept_id,
           year_of_birth
    FROM `{WORKSPACE_CDR}.person`
""").to_dataframe()

# GrafAnc sample IDs come from the FAM file — that's AoU person_id (string).
anc["person_id"] = anc["<sample_id_column e.g. IID>"].astype("int64")
merged = demo.merge(anc, on="person_id", how="inner")
merged.to_csv("grafanc_ancestry_with_demo.csv", index=False)

# Persist to the bucket.
os.system(
    f"gsutil cp grafanc_ancestry_with_demo.csv "
    f"{WORKSPACE_BUCKET}/grafanc/results/grafanc_ancestry_with_demo.csv"
)
