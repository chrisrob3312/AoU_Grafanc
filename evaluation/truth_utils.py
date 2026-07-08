"""
Shared helpers for reading haptools ground-truth local ancestry (.bp) and
evaluating a step function of ancestry along a chromosome. Imported by the
Tier-1 scorer (02) and the phenotype simulator (03) so the parsing lives in
one place (module name is a valid identifier, unlike the numbered scripts).
"""

import numpy as np

# Ancestry index -> label. MUST match the FLARE reference-panel label order.
ANC = {0: "EUR", 1: "AFR", 2: "AMR"}


def read_truth_bp(path):
    """
    Parse a haptools .bp file into {(sample, hap): [(end_bp, anc_index), ...]}
    sorted by position.

    haptools writes, per sample-haplotype block, a `Sample_<id>_<hap>` header
    then one line per segment: `pop  chrom  end_cM  end_bp`.
    """
    label_to_idx = {v: k for k, v in ANC.items()}
    truth, sample, hap = {}, None, 0
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("Sample"):
                token = line.split()[-1]
                *sid, h = token.rsplit("_", 1)
                sample, hap = "_".join(sid), int(h) - 1
                truth[(sample, hap)] = []
                continue
            pop, _chrom, _end_cm, end_bp = line.split()[:4]
            truth[(sample, hap)].append((int(end_bp), label_to_idx.get(pop, -1)))
    return truth


def truth_at_positions(segments, positions):
    """Ancestry of the segment whose end_bp first exceeds each position."""
    ends = np.array([e for e, _ in segments])
    ancs = np.array([a for _, a in segments])
    idx = np.clip(np.searchsorted(ends, positions, side="left"), 0, len(ancs) - 1)
    return ancs[idx]
