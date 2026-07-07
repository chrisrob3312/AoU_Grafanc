# Pipeline design notes & open decisions

Rationale, the choices baked into the scripts, and the things you should confirm
or decide before a production run.

## 1. Running GrafAnc on imputed + WGS (step 2–3)

- We take the **union** of participants across the imputed-array callset and the
  WGS callset, preferring WGS genotypes where a participant has both. This
  maximizes N (array-only participants are the majority of AoU) while keeping the
  highest-quality genotypes where available.
- GrafAnc needs the 282,424 ancestry SNPs. Imputed array data covers the vast
  majority of them at high quality; WGS covers all. GrafAnc still returns
  unbiased calls with some SNPs missing (≥10k for reliable continental, ≥50k for
  subcontinental), so mixing sources is fine — step 5 keeps the `#SNPs` count per
  participant and flags low-coverage ones.

### QC choices (step 3)

Deliberately light, per the GrafAnc docs (tolerates missingness, no LD pruning):

- **Kept:** biallelic SNPs, per-variant call rate ≥ 0.95, polymorphic; WGS
  entries gated on GQ ≥ 20 & DP ≥ 10; imputed sites assumed screened at R² ≥ 0.80
  upstream.
- **Not applied:** global HWE filter. The cohort is genetically structured, so
  HWE departures reflect ancestry, not error. Apply HWE only within a
  homogeneous subgroup if you want it at all.
- **Decision for you:** whether to also enforce imputation R² at extraction time
  for array-only participants (set `QC_IMPUTED_MIN_R2`, wire it into step 2's
  imputed read if the callset carries an INFO/R² field).

## 2. Cohort definitions (steps 5–7)

- Base cohorts come straight from `AncGroupID`. "LA1"/"LA2" in your notes = the
  Latin American 1/2 groups (601/602).
- **Multiracial (800) is the subtle one.** GrafAnc models only 3 continental
  components (European, African, East-Asian), so AMR admixture isn't directly
  readable. We resolve 800 with **supervised ADMIXTURE** over 7 superpops. This
  is why step 6 exists rather than reading AMR straight from GrafAnc.
- **Supervised mode also fixes cluster identity.** Running ADMIXTURE with the
  reference panels included and fixed labels (the `.pop` file) means each Q
  column *is* a named superpop — no post-hoc guessing which cluster is AMR vs
  AFR. Step 7 addresses columns by name via `superpop_column_order.txt`.
- **Thresholds** (step 7) are placeholders — set them to match your manuscript's
  admixture definitions. Current defaults:
  - 3-way: AFR+EUR+AMR ≥ 0.90, each ≥ 0.05;
  - 2-way: AFR+EUR ≥ 0.90, AMR < 0.05;
  - **plus an explicit cap on the "other" superpops** (`OTHER_MAX` = summed
    SAS+EAS+MEN+OCN ≤ 0.10, `OTHER_EACH_MAX` = no single one > 0.05). This is the
    "fourth/other component must be below X" gate you asked for — a participant
    with meaningful SAS/EAS/MEN/OCN ancestry is *not* a clean AFR-EUR(-AMR)
    mixture and is dropped to `multiracial_other`.
- A participant qualifying for both modes defaults to **3-way** (more general).
- **Marker set for step 6 is a flag** (`ancsnp` | `genomewide` | `both`). The
  ~282k ancestry SNPs are fine for global proportions and need no LD pruning
  (already selected for near-linkage-equilibrium). The `genomewide` path uses a
  denser, **LD-pruned** common-variant set (ADMIXTURE, unlike GrafAnc, assumes
  markers in linkage equilibrium) for tighter separation of the "other"
  superpops — which makes the `OTHER_MAX` gate more reliable. `both` runs each
  and feeds the genome-wide Q to step 7 by default (`PRIMARY_MARKERSET` overrides).

## 3. Reference panels (step 8)

- **Targeted panel:** EUR/AFR/AMR samples that are > 95% homogeneous per your
  local ADMIXTURE run, MXB included so the AMR reference has a genuine
  high-Amerindigenous anchor (typical AMR references are heavily admixed).
- **Comparison panel:** same EUR/AFR, broader/admixed AMR **without** MXB — the
  deliberately-less-fitting benchmark.
- FLARE reads whichever panel labels are in the ref map, so 2-way vs 3-way is
  just a matter of subsetting the map (EUR/AFR vs EUR/AFR/AMR) — step 10 does
  this automatically.
- **Pending:** the homogeneous sample lists from your HPC ADMIXTURE run. Until
  they land, the panel build can't run; everything upstream (GrafAnc,
  categorization) can.

## 4. FLARE (steps 9–10)

- Target genotypes for LAI are **genome-wide** markers (ACAF-threshold WGS /
  imputed), **not** the 282k ancestry SNPs — those were only for GrafAnc.
- Both target and reference must be **phased**. If AoU v9 ships phased WGS, skip
  SHAPEIT5 and point FLARE at those. Confirm.
- We run `probs=true` so the panel comparison (step 13) can use posterior
  sharpness.
- **Run the 2-way and 3-way individuals together (recommended).** `10_run_flare.sh
  combined <panel>` paints the union cohort against the 3-way EUR/AFR/AMR panel.
  FLARE assigns per-haplotype ancestry independently, so a genuinely 2-way
  AFR-EUR person just receives ~0 AMR — which is the correct answer — and you get
  one run and one output with no cohort split. Phase the combined list once:
  `09_phase_target.sh flare_admixed_all.samples.txt target_all`. The dedicated
  `2way` run (EUR/AFR-only panel) is kept only for the panel-fit comparison and
  for anyone who wants a strict no-AMR model.
- For the panel benchmark, each mode is still run against **both** the targeted
  and comparison panels so the comparison is apples-to-apples.

## 5. Post-processing (steps 11–13)

- **Uncalled variants (step 11):** local ancestry is piecewise-constant between
  recombination breakpoints, so an uncalled variant correctly inherits the
  ancestry of the tract it falls in. We collapse FLARE's per-marker calls into
  tracts, extend tract boundaries to inter-marker midpoints, and assign each
  target variant its enclosing tract — this is the "label by nearest tract" step
  you described.
- **Tractor (step 12):** emits the standard per-ancestry dosage + hapcount
  matrices Tractor's regression expects, with the FLARE tracts standing in for
  RFMix's MSP file.
- **Lab templates:** the two functions in step 11 and the one in step 12 are the
  designated plug-in points. The defaults are correct and runnable so the
  pipeline works today; swap in the lab versions if their tract definition or
  Tractor formatting differs.

## Open items / confirm before production

- [ ] v9 genomic env-var names & paths (WGS VDS, imputed MT, ACAF VCF).
- [ ] FLARE FORMAT field names (`AN1`/`AN2`/`ANP`) for your FLARE build.
- [ ] Whether AoU v9 provides pre-phased WGS (skip SHAPEIT5 if so).
- [ ] Admixture thresholds for multiracial → 2-way/3-way assignment.
- [ ] ADMIXTURE homogeneous sample lists (from HPC) + 7-superpop ref `.pop`.
- [ ] Lab tract-extraction / Tractor-prep templates.
