# Reference panel — homogeneous samples by super-population

The LAI reference panel is the **homogeneous** per-super-population samples from
the local ADMIXTURE homogeneity screen (`04_homogeneity_panel`), NOT the whole
callset. Five super-populations:

| Super-pop | N | Source |
|-----------|---|--------|
| AFR | 634 | HGDP + 1KG |
| AMR | 88  | HGDP + 1KG + **MX Biobank** (see below) |
| EAS | 667 | HGDP + 1KG |
| EUR | 620 | HGDP + 1KG |
| SAS | 48  | HGDP |
| **total** | **2057** | |

Files: `homog_by_super_pop/{AFR,AMR,EAS,EUR,SAS}_ids.txt`, one ID per line.

These are **de-identified sample IDs / metadata only — no genetic sequences.**
They are committed so that anyone with the relevant data access (1KG/HGDP are
public; MX Biobank via its DUA) can reproduce the exact panel composition.

## The MXB contribution is built into the AMR list

The **AMR list of 88 = 50 MX Biobank (`MXB_*`) + 38 HGDP/1KG Amerindigenous.**
Dropping the `MXB_*` IDs isolates the MXB contribution exactly, which defines the
two headline panel arms:

| Arm | AMR reference | N (AFR+EUR+AMR) | Isolates |
|-----|---------------|-----------------|----------|
| `mxb_expanded` | all 88 (incl. 50 MXB) | 1342 | the proposal — large **and** matched |
| `amr_small_homog` | 38 (MXB removed) | 1292 | matching alone, low N |

The other two evaluation arms use AMR references **outside** this homogeneity
panel and must be supplied separately:
- `amr_large_admixed` — full admixed 1KG AMR (MXL/PEL/CLM/PUR): size, contaminated;
- `aou_default_1kg` — 1KG-style baseline ≈ what AoU CDR v9 will likely ship.

## Build

```bash
bash reference/build_panels.sh
```
Produces (git-ignored, regenerate anytime) in `reference/panels/`:
- `ref_labels.tsv` — `sample_id <TAB> superpop` for RFMix classes / FLARE ref map / supervised ADMIXTURE
- `AMR_mxb.txt`, `AMR_nonmxb.txt` — the 50/38 split
- `mxb_expanded.keep`, `amr_small_homog.keep` — LAI panel arms

Copy `reference/homog_by_super_pop/` (or the built keep-lists) into the workspace
bucket for the panel-build step (`flare/ref_panel/08_build_ref_panels.sh`).
