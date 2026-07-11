"""
Shared helpers for reading admix-simu ground-truth local ancestry and evaluating
a step function of ancestry along a chromosome. Imported by the Tier-1 scorer
(02) and the phenotype simulator (03).

admix-simu (`simu-mix.pl`) writes a `.bp` breakpoint file: two lines per
simulated individual (haplotype 1 then haplotype 2), each line a space-separated
list of `POP:END` tokens giving the ancestry and END coordinate of each segment
along the chromosome, in order. Sample IDs come from a companion order file
(one ID per individual), because the `.bp` itself is positional.

>>> CONFIRM against your admix-simu fork <<<
Two things vary between forks and MUST match your build:
  * END is physical bp here (set END_UNIT="cM" + pass a map if yours is genetic).
  * Two consecutive lines per individual = hap1, hap2 (set HAPS_INTERLEAVED).
"""

import numpy as np

# Ancestry index -> label. MUST match the LAI reference-panel label order.
ANC = {0: "EUR", 1: "AFR", 2: "AMR"}
_LABEL_TO_IDX = {v: k for k, v in ANC.items()}

END_UNIT = "bp"          # "bp" or "cM"
HAPS_INTERLEAVED = True   # True: line 2i, 2i+1 are the two haps of individual i


def read_truth_bp(bp_path, ids_path):
    """
    Parse admix-simu `.bp` + a sample-order file into
    {(sample, hap): [(end_bp, anc_index), ...]} sorted by position.
    """
    sample_ids = [l.strip() for l in open(ids_path) if l.strip()]
    lines = [l.rstrip("\n") for l in open(bp_path) if l.strip()]
    if HAPS_INTERLEAVED and len(lines) != 2 * len(sample_ids):
        raise ValueError(
            f"{len(lines)} .bp lines vs {len(sample_ids)} ids "
            f"(expected 2x for interleaved haplotypes) — check the format flags."
        )

    truth = {}
    for i, sample in enumerate(sample_ids):
        for hap in (0, 1):
            line = lines[2 * i + hap] if HAPS_INTERLEAVED else lines[hap * len(sample_ids) + i]
            segs = []
            for tok in line.split():
                pop, end = tok.split(":")
                segs.append((int(float(end)), _LABEL_TO_IDX.get(pop, -1)))
            truth[(sample, hap)] = sorted(segs)
    return truth


def truth_at_positions(segments, positions):
    """Ancestry of the segment whose END first reaches each position."""
    ends = np.array([e for e, _ in segments])
    ancs = np.array([a for _, a in segments])
    idx = np.clip(np.searchsorted(ends, positions, side="left"), 0, len(ancs) - 1)
    return ancs[idx]
