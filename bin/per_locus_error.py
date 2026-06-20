#!/usr/bin/env python3
"""Per-locus sequencing-error rates from a samtools mpileup.

Computes, for each candidate locus, the mean Phred-scaled base-error probability
across all bases piled up at that site (all individuals pooled). Output is one
value per locus, one per line, in the order given by ``--loci`` so it lines up
with the EBG read-count matrices (which are built from the same VCF).

Modelled on EBG's helper ``per-locus-err.py`` but driven by the candidate-loci
list rather than pileup line order, so loci with no coverage (absent from the
pileup) still get a value (the default) and alignment is never broken.

mpileup column layout: chrom, pos, ref, then per individual a triple of
(depth, read-bases, base-qualities). The quality strings are the 6th column
onward, every third field.

Usage: per_locus_error.py --pileup <file> --loci <chrom\\tpos file> [--default 0.01]
"""

import argparse
import sys


def phred_error(qual_char):
    return 10.0 ** (-(ord(qual_char) - 33) / 10.0)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pileup", required=True, help="samtools mpileup output")
    parser.add_argument("--loci", required=True, help="CHROM<TAB>POS per locus, in matrix order")
    parser.add_argument("--default", type=float, default=0.01, help="Error used for loci with no quality data")
    parser.add_argument("--out", default="error.txt", help="Output file (default error.txt)")
    args = parser.parse_args(argv)

    # locus key "chrom:pos" -> mean error
    err_by_locus = {}
    with open(args.pileup) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            chrom, pos = fields[0], fields[1]
            quals = "".join(fields[5::3])  # 6th col onward, every 3rd = qualities
            if quals:
                errs = [phred_error(c) for c in quals]
                err_by_locus[f"{chrom}:{pos}"] = sum(errs) / len(errs)

    with open(args.loci) as loci_fh, open(args.out, "w") as out_fh:
        for line in loci_fh:
            line = line.strip()
            if not line:
                continue
            chrom, pos = line.split("\t")[:2]
            out_fh.write(f"{err_by_locus.get(f'{chrom}:{pos}', args.default):.6g}\n")

    return 0


if __name__ == "__main__":
    sys.exit(main())
