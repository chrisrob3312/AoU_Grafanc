# GRAF-anc on the All of Us Researcher Workbench (v8)

Pipeline to assign global ancestry to every AoU v8 participant with WGS data
using the ~280 ancestry-informative SNPs (AIMs) from GrafAnc
(<https://github.com/jimmy-penn/GrafAnc>  ← replace with actual URL).

## Workflow

| Step | File | Where to run |
|------|------|--------------|
| 1. Pull GrafAnc and extract the AIM SNP list | `scripts/01_get_grafanc_snps.sh` | RW Cloud Analysis terminal |
| 2. Subset AoU WGS VDS to the AIM SNPs (Hail) | `notebooks/02_extract_aou_snps.py` | Dataproc / Hail env |
| 3. Intersect with the lab's v8 QC-pass list, export PLINK | `notebooks/03_filter_qc_and_export_plink.py` | Dataproc / Hail env |
| 4. Run GrafAnc on the PLINK fileset | `scripts/04_run_grafanc.sh` | RW Cloud Analysis terminal |
| 5. Parse + summarize ancestry calls | `notebooks/05_parse_results.py` | Standard Jupyter env |

## Things you fill in

Anywhere you see `<...>` you need to drop in:
- `<copy and paste workspace name>` – the AoU workspace this runs in
- `<path/to/labmate/qc_pass_variants.tsv>` – the QC-pass variant list
- `<grafanc_repo_url>` – Jimmy's GitHub URL
- `<num_workers>` – Dataproc worker count

All other AoU paths are read from the workbench environment variables
(`WORKSPACE_BUCKET`, `WGS_VDS_PATH`, `WORKSPACE_CDR`, ...) so the code is
portable across workspaces.
