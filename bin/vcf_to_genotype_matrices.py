#!/usr/bin/env python3
"""Build EBG and updog read-count matrices for one ploidy group from a VCF.

Reads a (bgzipped or plain) biallelic VCF that carries per-genotype ``AD`` and
writes, for the subset of samples whose ploidy matches ``--ploidy``:

  <prefix>.tot.txt    EBG total reads   (rows = individuals, cols = loci; space
                      separated; missing = -9)
  <prefix>.alt.txt    EBG ALT reads     (same layout as tot)
  <prefix>.ref.txt    updog ref counts  (rows = loci, cols = individuals; TSV
                      with a header row of sample ids and a first column of SNP
                      ids; missing = 0)
  <prefix>.size.txt   updog total reads (same layout as ref; missing = 0)
  <prefix>.samples.txt  sample ids, one per line, in EBG row / updog column order
  <prefix>.loci.txt     "CHROM\\tPOS" per locus, in EBG column / updog row order
  <prefix>.nind         the number of individuals in this group

EBG encodes missing genotypes as -9; updog treats missing as zero reads (size 0).
Loci are emitted in VCF (coordinate) order so the EBG matrices line up with the
per-locus error vector, which is built from the same VCF.
"""

import argparse
import gzip
import sys


def open_maybe_gz(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def load_ploidy_map(path):
    """sample -> ploidy (int) from a 'sample<TAB>ploidy' file."""
    mapping = {}
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            sample, ploidy = line.split("\t")[:2]
            mapping[sample] = int(ploidy)
    return mapping


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", required=True, help="Candidate VCF (.vcf or .vcf.gz)")
    parser.add_argument("--ploidy-map", required=True, help="sample<TAB>ploidy file")
    parser.add_argument("--ploidy", type=int, required=True, help="Target ploidy for this group")
    parser.add_argument("--prefix", required=True, help="Output file prefix")
    args = parser.parse_args(argv)

    ploidy_of = load_ploidy_map(args.ploidy_map)

    cols = None        # indices (into the sample columns) kept for this group
    samples = None     # sample ids kept for this group
    loci = []          # (chrom, pos)
    tot = []           # tot[locus][sample_in_group]
    alt = []

    with open_maybe_gz(args.vcf) as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                all_samples = fields[9:]
                cols = [
                    i for i, s in enumerate(all_samples)
                    if ploidy_of.get(s) == args.ploidy
                ]
                samples = [all_samples[i] for i in cols]
                continue
            if cols is None:
                sys.exit("error: VCF had no #CHROM header line")
            fmt = fields[8].split(":")
            try:
                ad_idx = fmt.index("AD")
            except ValueError:
                sys.exit("error: FORMAT has no AD field; cannot extract read counts")

            loci.append((fields[0], fields[1]))
            tot_row, alt_row = [], []
            for i in cols:
                gt = fields[9 + i]
                subfields = gt.split(":")
                ad = subfields[ad_idx] if len(subfields) > ad_idx else "."
                if gt.startswith(".") or ad in (".", "") or "." in ad.split(","):
                    tot_row.append(-9)
                    alt_row.append(-9)
                else:
                    parts = ad.split(",")
                    ref_c = int(parts[0])
                    alt_c = int(parts[1]) if len(parts) > 1 else 0
                    tot_row.append(ref_c + alt_c)
                    alt_row.append(alt_c)
            tot.append(tot_row)
            alt.append(alt_row)

    if not samples:
        sys.exit(f"error: no samples with ploidy {args.ploidy} found in the VCF")

    n_loci = len(loci)
    n_ind = len(samples)

    # EBG matrices: rows = individuals, cols = loci, space separated, missing -9.
    with open(f"{args.prefix}.tot.txt", "w") as ft, open(f"{args.prefix}.alt.txt", "w") as fa:
        for j in range(n_ind):
            ft.write(" ".join(str(tot[l][j]) for l in range(n_loci)) + "\n")
            fa.write(" ".join(str(alt[l][j]) for l in range(n_loci)) + "\n")

    # updog matrices: rows = loci, cols = individuals; TSV with names; missing -> 0.
    snp_ids = [f"{c}_{p}" for c, p in loci]
    header = "snp\t" + "\t".join(samples) + "\n"
    with open(f"{args.prefix}.ref.txt", "w") as fr, open(f"{args.prefix}.size.txt", "w") as fs:
        fr.write(header)
        fs.write(header)
        for l in range(n_loci):
            ref_cells, size_cells = [], []
            for j in range(n_ind):
                t = tot[l][j]
                a = alt[l][j]
                if t < 0:
                    ref_cells.append("0")
                    size_cells.append("0")
                else:
                    ref_cells.append(str(t - a))
                    size_cells.append(str(t))
            fr.write(snp_ids[l] + "\t" + "\t".join(ref_cells) + "\n")
            fs.write(snp_ids[l] + "\t" + "\t".join(size_cells) + "\n")

    with open(f"{args.prefix}.samples.txt", "w") as f:
        f.write("\n".join(samples) + "\n")
    with open(f"{args.prefix}.loci.txt", "w") as f:
        f.write("\n".join(f"{c}\t{p}" for c, p in loci) + "\n")
    with open(f"{args.prefix}.nind", "w") as f:
        f.write(f"{n_ind}\n")

    sys.stderr.write(
        f"[vcf_to_genotype_matrices] ploidy {args.ploidy}: "
        f"{n_ind} individuals x {n_loci} loci\n"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
