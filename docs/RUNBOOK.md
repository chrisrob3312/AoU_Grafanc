# Runbook — what to provide, and how to run inside the AoU Workbench

This is the operational companion to `README.md` (what the pipeline is) and
`docs/PIPELINE.md` (why each choice was made). It covers:

1. **Pre-flight checklist** — everything you must hand over / produce first.
2. **Getting your data and the tools into the Verily Workbench.**
3. **Run order**, with the environment each step needs.

> You code in shell/R; several steps are Python (Hail has no R API, and FLARE
> post-processing is easiest in Python). Every Python file has a plain-language
> header explaining what it does and heavy inline comments. The pure-tabular
> steps (5, 7, 13) *could* be ported to R if you'd prefer — say the word.

---

## 1. Pre-flight checklist

Nothing runs end-to-end until these are in place. Grouped by where they come from.

### A. From your local HPC ADMIXTURE run (pending)
- [ ] `eur_homog95.samples.txt` — EUR reference samples >95% homogeneous (one ID/line)
- [ ] `afr_homog95.samples.txt` — AFR reference samples >95% homogeneous
- [ ] `amr_homog95.samples.txt` — AMR reference samples >95% homogeneous (incl. MXB)
- [ ] `amr_broad.samples.txt` — broader/admixed AMR set **without** MXB (comparison panel)

### B. Reference data to upload into the workspace bucket
- [ ] `1kg_hgdp_mxb.joint.vcf.gz` (+ `.tbi`) — your joint-called combined callset
- [ ] 7-superpop reference PLINK sets for step 6:
      - `7superpop_ref.ancsnp.{bed,bim,fam}` (restricted to the ancestry SNPs)
      - `7superpop_ref.genomewide.{bed,bim,fam}` (only if you run the `genomewide` marker set)
- [ ] `7superpop_ref_labels.tsv` — two columns: `sample_id <TAB> superpop`
      (superpop ∈ AMR/AFR/EUR/SAS/EAS/MEN/OCN), reference samples only
- [ ] GRCh38 genetic maps for phasing + FLARE (`GENETIC_MAP_DIR`)

### C. Lab templates (pending — plug-in points already stubbed)
- [ ] FLARE tract-extraction template → `flare/postprocess/11_fill_uncalled_tracts.py`
      (`flare_to_tracts()` / `assign_variant_ancestry()`)
- [ ] Tractor-prep template → `flare/postprocess/12_flare_to_tractor.py`
      (`write_tractor_matrices()`)

### D. Config values to confirm in `config/config.sh`
- [ ] `WORKSPACE_NAME`
- [ ] v9 genomic env-var names/paths: `WGS_VDS_PATH`, `IMPUTED_MT_PATH`, and the
      per-chrom ACAF VCF used in steps 6/9 (`ACAF_VCF`). **Confirm these against
      the "Genomic data" page in your v9 workbench** — v8 names are the fallback.
- [ ] FLARE FORMAT field names (`AN1`/`AN2`/`ANP`) — check the header of the
      `flare.jar` you download; they're flagged inline in steps 11–13.
- [ ] Thresholds in step 7 (`THREE_WAY_MIN_SUM`, `EACH_MIN`, `TWO_WAY_MIN_SUM`,
      `AMR_MAX`, `OTHER_MAX`, `OTHER_EACH_MAX`) — set to your manuscript's defs.
- [ ] Compute sizing: `N_THREADS`, `N_WORKERS` (Dataproc), FLARE `-Xmx` heap,
      GrafAnc `--maxmem`.

### E. Decisions to make
- [ ] Step 6 marker set: `ancsnp` (fast), `genomewide` (tighter "other"
      separation), or `both`.
- [ ] FLARE cohort mode: `combined` (recommended — everyone vs the 3-way panel)
      vs separate `2way`/`3way` runs.
- [ ] Whether AoU v9 ships **pre-phased** WGS (if so, skip SHAPEIT5 in steps 8–9).

---

## 2. Getting data and tools into the Workbench

### 2.1 How AoU data is already available (don't upload these)
- **Phenotypes / CDR** live in BigQuery; the dataset is in `$WORKSPACE_CDR`.
- **Genomic data** (WGS VDS, imputed callset, ACAF VCFs) live in Google Cloud
  Storage and are addressed through environment variables AoU sets in every
  analysis environment (e.g. `$WGS_VDS_PATH`). You read them in place — no copy.
- Confirm the exact v9 variable names under **"Genomic data" → "How to use"** in
  the workbench; drop them into `config/config.sh`.

### 2.2 Uploading YOUR reference panels + sample lists (checklist A & B)
Your combined 1KG-HGDP+MXB callset, the 7-superpop reference sets, the ADMIXTURE
sample lists, and the genetic maps are your own files and go into the
**workspace bucket** (`$WORKSPACE_BUCKET`), which persists across sessions.

Two ways in:

**(a) Straight from your HPC to the bucket (best for the big VCFs).**
On a machine that has your files and the gcloud SDK + workbench credentials:
```bash
# One-time: authenticate to the account tied to your AoU workbench.
gcloud auth login
# Copy (multi-threaded) into the workspace bucket. Get the bucket path by
# running:  echo $WORKSPACE_BUCKET   inside a workbench terminal.
gsutil -m cp 1kg_hgdp_mxb.joint.vcf.gz* \
    gs://<your-workspace-bucket>/grafanc_flare_v9/flare/ref_panels/combined/
gsutil -m cp *_homog95.samples.txt amr_broad.samples.txt \
    gs://<your-workspace-bucket>/grafanc_flare_v9/cohorts/
gsutil -m cp 7superpop_ref.*.{bed,bim,fam} 7superpop_ref_labels.tsv \
    gs://<your-workspace-bucket>/grafanc_flare_v9/ref/
```

**(b) Through the Jupyter UI (fine for the small text files).**
In the workbench, open Jupyter → **Upload** → then move the file into the bucket
from a terminal with `gsutil cp <file> $WORKSPACE_BUCKET/...`. The local disk is
ephemeral; only the bucket survives, so always land files in `$WORKSPACE_BUCKET`.

> **Data-use note:** bringing your own MXB reference data into the workspace must
> comply with the AoU Data User Code of Conduct and your MXB data agreement.
> Reference panels stay in your private workspace bucket. This pipeline never
> exports AoU participant-level data out of the workbench.

### 2.3 Installing the tools
Do this in a **Cloud Analysis** environment (a "General Analysis" /
Jupyter+terminal VM). Install into `~/tools` and stage anything reusable in the
bucket. Reinstalls are cheap; the scripts assume these are on `PATH`.

| Tool | How | Used by |
|------|-----|---------|
| **GrafAnc** | `git clone` + `make` — done for you by `scripts/01_setup_grafanc.sh` | steps 1, 4 |
| **PLINK 1.9 / 2** | `conda install -c bioconda plink plink2` (or download binaries) | steps 3, 6, 8 |
| **bcftools** | `conda install -c bioconda bcftools` | steps 6, 8, 9 |
| **ADMIXTURE** | download the static binary from the Alkes/UCLA page, `chmod +x` | step 6 |
| **SHAPEIT5** | `conda install -c bioconda shapeit5` (or binaries) — skip if v9 WGS is pre-phased | steps 8, 9 |
| **FLARE** | `wget` the `flare.jar` (needs Java 8+: `conda install -c conda-forge openjdk`) — done by `flare/10_run_flare.sh` | step 10 |
| **Tractor** | `git clone https://github.com/Atkinson-Lab/Tractor` | after step 12 |
| **Hail** | pre-installed on the Dataproc/Hail environment | steps 2, 3 |

A quick one-shot install for the terminal tools:
```bash
conda install -y -c bioconda -c conda-forge \
    plink plink2 bcftools shapeit5 openjdk
```

### 2.4 Which environment runs which step
- **Cloud Analysis terminal** (CPU VM): steps 1, 4, 6, 8, 9, 10.
- **Dataproc / Hail** cluster: steps 2, 3. Spin up with `$N_WORKERS` workers.
- **Plain Jupyter (pandas)**: steps 5, 7, 11, 12, 13.

---

## 3. Run order

```bash
# --- one-time: edit config, then source it in every terminal session ---------
vim config/config.sh                 # fill in the <...> placeholders
source config/config.sh

# --- GrafAnc global ancestry (Cloud Analysis terminal + Hail) -----------------
./scripts/01_setup_grafanc.sh                 # clone/build GrafAnc, emit SNP list
#   run notebooks/02_extract_anc_snps.py       on Dataproc/Hail
#   run notebooks/03_qc_and_export.py          on Dataproc/Hail
./scripts/04_run_grafanc.sh                    # run GrafAnc → per-participant AncGroupID

# --- categorize + resolve multiracial ----------------------------------------
#   run notebooks/05_categorize_cohorts.py     (Jupyter)
./admixture/06_multiracial_supervised_admixture.sh both     # ancsnp|genomewide|both
#   run notebooks/07_define_flare_cohorts.py   (Jupyter) → final 2way/3way/all lists

# --- FLARE reference panels + target phasing ---------------------------------
./flare/ref_panel/08_build_ref_panels.sh                    # targeted + comparison panels
./flare/09_phase_target.sh flare_admixed_all.samples.txt target_all

# --- run FLARE (combined cohort vs each panel) -------------------------------
./flare/10_run_flare.sh combined targeted
./flare/10_run_flare.sh combined comparison
#   (optional strict runs for the benchmark: 2way targeted / 3way targeted, etc.)

# --- post-process + Tractor prep + panel comparison (Jupyter) -----------------
#   run flare/postprocess/11_fill_uncalled_tracts.py  (per chrom)
#   run flare/postprocess/12_flare_to_tractor.py      (per chrom)
#   run flare/compare/13_panel_accuracy_comparison.py
```

Steps 2, 3, 5, 7, 11, 12, 13 are notebooks — open each `.py` in Jupyter (or
`%run` it) rather than executing from the shell, so you can inspect the printed
counts at each stage. Every one echoes row/sample counts so you can sanity-check
before moving on.
