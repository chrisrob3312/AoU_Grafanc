"""
Step 11 — propagate FLARE local-ancestry calls to variants FLARE did not
directly output.

FLARE emits a per-marker local-ancestry call (per haplotype) at the markers it
modeled.  Downstream tools (and any variant-level analysis) need an ancestry
label at EVERY target variant.  Because local ancestry is piecewise-constant
along a chromosome (it changes only at recombination breakpoints / tract
boundaries), the correct fill is: each uncalled variant inherits the ancestry
of the tract it falls in — i.e. the nearest called marker on the same
haplotype, within the enclosing ancestry segment.

This script:
  1. reads FLARE's ANC output (per-haplotype ancestry codes at called markers),
  2. reconstructs continuous ancestry TRACTS per haplotype,
  3. assigns every target variant position the ancestry of its enclosing tract
     (nearest-marker tie-break at exact boundaries),
  4. writes a per-variant, per-haplotype ancestry table ready for step 12.

>>> LAB TEMPLATE PLUG-IN <<<
Your lab's "extract tracts from FLARE" and "post-processing" templates slot into
`flare_to_tracts()` and `assign_variant_ancestry()`.  The reference
implementation below is correct and runnable; swap in the lab versions when they
arrive if their tract definition differs.
"""

import argparse
import numpy as np
import pandas as pd
import pysam   # cyvcf2 also fine; pysam reads FLARE's .anc.vcf.gz


def flare_to_tracts(flare_anc_vcf):
    """
    Read FLARE's .anc.vcf.gz and return, per sample-haplotype, the ordered list
    of ancestry calls at each called marker:
        { (sample, hap): DataFrame[pos, anc] }  (hap in {0,1})

    FLARE encodes the per-haplotype ancestry in the FORMAT 'AN1'/'AN2' fields
    (integer ancestry index matching the ref-panel label order).  Confirm the
    field names against your FLARE version's header.
    """
    vcf = pysam.VariantFile(flare_anc_vcf)
    samples = list(vcf.header.samples)
    rows = {(s, h): [] for s in samples for h in (0, 1)}
    for rec in vcf:
        for s in samples:
            call = rec.samples[s]
            an1 = call.get("AN1"); an2 = call.get("AN2")   # <-- confirm field names
            rows[(s, 0)].append((rec.pos, an1))
            rows[(s, 1)].append((rec.pos, an2))
    return {k: pd.DataFrame(v, columns=["pos", "anc"]) for k, v in rows.items()}


def build_segments(marker_df):
    """
    Collapse consecutive equal-ancestry markers into continuous [start,end,anc]
    tracts.

    Plain terms: FLARE gives us an ancestry label at each marker along the
    chromosome, e.g.
        pos 100=AFR, 200=AFR, 300=EUR, 400=EUR, 500=AFR
    Local ancestry only changes at recombination breakpoints, so runs of the
    same label are one tract. We turn the per-marker list into:
        AFR from 100-250, EUR from 251-450, AFR from 451-end
    (the 250/450 cut points are midpoints between the last marker of one tract
    and the first of the next, so every base in between belongs to exactly one
    tract and there are no gaps).
    """
    # Drop markers with no call and sort by position (defensive; FLARE is sorted).
    m = marker_df.dropna().sort_values("pos").reset_index(drop=True)
    if m.empty:
        return pd.DataFrame(columns=["start", "end", "anc"])

    # `change` increments every time the ancestry label differs from the previous
    # row, giving each run of identical labels a unique group number. Grouping by
    # it collapses each run into one row: its min pos, max pos, and the label.
    change = m["anc"].ne(m["anc"].shift()).cumsum()
    segs = m.groupby(change).agg(start=("pos", "min"),
                                 end=("pos", "max"),
                                 anc=("anc", "first")).reset_index(drop=True)

    # Now widen each tract so consecutive tracts touch with no gap. For each
    # adjacent pair, the boundary is the midpoint between one tract's last marker
    # (`end`) and the next tract's first marker (`start`).
    mids = ((segs["end"].values[:-1] + segs["start"].values[1:]) / 2).astype(int)
    segs.loc[:len(mids) - 1, "end"] = mids        # each tract ends at the midpoint...
    segs.loc[1:, "start"] = mids + 1              # ...and the next starts just after
    segs.loc[0, "start"] = 0                       # first tract runs from the chrom start
    segs.loc[len(segs) - 1, "end"] = np.iinfo(np.int64).max  # last runs to the chrom end
    return segs


def assign_variant_ancestry(segments, variant_positions):
    """
    For each target variant position, return the ancestry of the tract it sits in.

    Plain terms: given the tract table from build_segments (sorted, gap-free) and
    a list of variant positions, find which tract each position falls in and
    return that tract's ancestry label. This is the "label the uncalled variant
    by the tract it's closest to / inside" step.
    """
    starts = segments["start"].values
    # np.searchsorted finds, for each variant position, where it would slot into
    # the sorted list of tract start positions. "-1" turns that into the index of
    # the tract whose start is <= the position, i.e. the enclosing tract.
    idx = np.searchsorted(starts, variant_positions, side="right") - 1
    idx = np.clip(idx, 0, len(segments) - 1)      # guard the edges
    return segments["anc"].values[idx]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--flare_anc_vcf", required=True, help="FLARE .anc.vcf.gz for one chrom")
    ap.add_argument("--target_sites", required=True,
                    help="TSV/PLINK bim with every target variant position for this chrom")
    ap.add_argument("--out", required=True, help="output per-variant ancestry table (.tsv.gz)")
    args = ap.parse_args()

    # Target variant positions (one chromosome).
    sites = pd.read_csv(args.target_sites, sep="\t")
    positions = sites["pos"].values if "pos" in sites.columns else sites.iloc[:, 3].values

    tracts_by_hap = flare_to_tracts(args.flare_anc_vcf)

    out = {"pos": positions}
    for (sample, hap), markers in tracts_by_hap.items():
        segs = build_segments(markers)
        if segs.empty:
            out[f"{sample}_hap{hap}"] = np.nan
        else:
            out[f"{sample}_hap{hap}"] = assign_variant_ancestry(segs, positions)

    pd.DataFrame(out).to_csv(args.out, sep="\t", index=False, compression="gzip")
    print(f"Wrote per-variant ancestry: {args.out}  "
          f"({len(positions):,} variants × {len(tracts_by_hap)} haplotypes)")


if __name__ == "__main__":
    main()
