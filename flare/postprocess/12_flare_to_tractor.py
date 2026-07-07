"""
Step 12 — convert filled per-variant local-ancestry calls into Tractor inputs.

Tractor (https://github.com/Atkinson-Lab/Tractor) runs ancestry-aware GWAS.  For
K ancestries it needs, at every variant, per individual:
    * an ancestry-DOSAGE count  (0/1/2 copies from each ancestry), and
    * ancestry-specific haplotype counts / genotypes.

Tractor's `ExtractTracts.py` normally derives these from a phased VCF + the MSP
tracts from RFMix.  Here the FLARE-derived per-haplotype ancestry table
(step 11) plays the role of the MSP file, so we emit the same
`*.ancK.dosage.txt` / `*.hapcountK.txt` matrices Tractor's regression expects.

>>> LAB TEMPLATE PLUG-IN <<<
Drop your lab's "get FLARE output ready for Tractor" template into
`write_tractor_matrices()`.  The implementation below produces the standard
Tractor per-ancestry dosage + hapcount files and is a valid default.
"""

import argparse
import numpy as np
import pandas as pd
import pysam

# Ancestry index -> label MUST match the FLARE ref-panel label order.
ANCESTRY_LABELS = {0: "EUR", 1: "AFR", 2: "AMR"}   # <-- confirm vs your ref map


def load_phased_haplotypes(phased_vcf, chrom):
    """Return a DataFrame [pos] + one column per sample_hap of {0,1} allele calls."""
    vcf = pysam.VariantFile(phased_vcf)
    samples = list(vcf.header.samples)
    data = {"pos": []}
    for s in samples:
        data[f"{s}_hap0"] = []
        data[f"{s}_hap1"] = []
    for rec in vcf:
        data["pos"].append(rec.pos)
        for s in samples:
            a = rec.samples[s].allele_indices    # phased (a0, a1)
            data[f"{s}_hap0"].append(a[0])
            data[f"{s}_hap1"].append(a[1])
    return pd.DataFrame(data)


def write_tractor_matrices(hap_alleles, hap_ancestry, out_prefix, labels=ANCESTRY_LABELS):
    """
    For each ancestry k, write:
      {out_prefix}.anc{k}.dosage.txt  — per-individual copies of ALT from anc k
      {out_prefix}.hapcount{k}.txt    — per-individual # haplotypes from anc k
    aligned on the intersection of variant positions.
    """
    merged = hap_alleles.merge(hap_ancestry, on="pos", suffixes=("_al", "_an"))
    pos = merged["pos"].values
    samples = sorted({c.rsplit("_hap", 1)[0] for c in hap_alleles.columns if c != "pos"})

    # Plain terms: for each ancestry k (EUR/AFR/AMR) we build two matrices, one
    # row per variant, one column per individual:
    #   hapcount = how many of the person's 2 haplotypes are of ancestry k here
    #              (0, 1, or 2)
    #   dosage   = how many ALT alleles the person carries that sit on an
    #              ancestry-k haplotype here (0, 1, or 2)
    # Tractor's regression uses exactly these per-ancestry counts.
    for k, lab in labels.items():
        dosage = {"pos": pos}
        hapcount = {"pos": pos}
        for s in samples:
            dos = np.zeros(len(pos)); hc = np.zeros(len(pos))
            for hap in (0, 1):
                # allele: 0=REF, 1=ALT on this haplotype at each variant.
                allele = merged[f"{s}_hap{hap}_al"].values
                # anc: which ancestry index this haplotype was painted at each variant.
                anc    = merged[f"{s}_hap{hap}_an"].values
                from_k = (anc == k)                            # True where this hap is ancestry k
                hc  += from_k.astype(int)                      # count the haplotype
                dos += ((allele == 1) & from_k).astype(int)    # count ALT copies on ancestry-k haps
            dosage[s] = dos.astype(int)
            hapcount[s] = hc.astype(int)
        pd.DataFrame(dosage).to_csv(f"{out_prefix}.anc{k}_{lab}.dosage.txt",
                                    sep="\t", index=False)
        pd.DataFrame(hapcount).to_csv(f"{out_prefix}.hapcount{k}_{lab}.txt",
                                      sep="\t", index=False)
        print(f"  wrote {lab}: dosage + hapcount ({len(pos):,} variants)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--phased_vcf", required=True, help="phased target VCF for this chrom")
    ap.add_argument("--ancestry_table", required=True,
                    help="per-variant per-hap ancestry table from step 11 (.tsv.gz)")
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--out_prefix", required=True)
    args = ap.parse_args()

    hap_alleles = load_phased_haplotypes(args.phased_vcf, args.chrom)
    hap_ancestry = pd.read_csv(args.ancestry_table, sep="\t")
    write_tractor_matrices(hap_alleles, hap_ancestry, args.out_prefix)
    print(f"Tractor inputs → {args.out_prefix}.*")


if __name__ == "__main__":
    main()
