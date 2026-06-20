#!/usr/bin/env python3
"""Compare EBG and updog polyploid genotype calls and report concordance.

Both tools are run once per ploidy group; this script discovers each group from
the files in the working directory, harmonises the two dosage conventions, and
writes a per-call table plus a concordance summary.

Per group ``g`` it expects (all produced upstream, named by ploidy):
  ebg_p<g>-genos.txt   EBG calls, rows = individuals, cols = loci, ALT-allele
                       dosage 0..g (missing = -9)
  p<g>.samples.txt     sample ids in EBG row order
  p<g>.loci.txt        "CHROM\\tPOS" in EBG column order
  updog_p<g>.updog.tsv  long table: snp, ind, geno (REFERENCE-allele dosage)

EBG reports ALT dosage; updog reports REFERENCE dosage, so updog is converted to
ALT dosage as ``ploidy - geno`` before comparison. SNP ids are ``CHROM_POS`` on
both sides.

Outputs:
  genotype_concordance.tsv          snp, ind, ploidy, ebg_alt, updog_alt, status
  genotype_concordance_summary.tsv  overall / per-ploidy / per-dosage agreement
"""

import glob
import os
import re
import sys
from collections import defaultdict


def load_samples(path):
    with open(path) as fh:
        return [line.strip() for line in fh if line.strip()]


def load_loci(path):
    snps = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line:
                chrom, pos = line.split("\t")[:2]
                snps.append(f"{chrom}_{pos}")
    return snps


def load_ebg(path, samples, snps):
    """Return {(snp, ind): alt_dosage or None for missing}."""
    calls = {}
    with open(path) as fh:
        rows = [ln.rstrip("\n").rstrip("\t").split("\t") for ln in fh if ln.strip()]
    if len(rows) != len(samples):
        sys.exit(f"error: {path} has {len(rows)} rows but {len(samples)} samples")
    for i, row in enumerate(rows):
        if len(row) != len(snps):
            sys.exit(f"error: {path} row {i} has {len(row)} cols but {len(snps)} loci")
        for j, val in enumerate(row):
            v = int(val)
            calls[(snps[j], samples[i])] = None if v < 0 else v
    return calls


def load_updog(path, ploidy):
    """Return {(snp, ind): alt_dosage or None}, converting ref->alt dosage."""
    calls = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {name: k for k, name in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            snp, ind, geno = f[idx["snp"]], f[idx["ind"]], f[idx["geno"]]
            if geno in ("NA", ""):
                calls[(snp, ind)] = None
            else:
                calls[(snp, ind)] = ploidy - int(round(float(geno)))
    return calls


def main():
    groups = sorted(
        int(re.search(r"ebg_p(\d+)-genos\.txt$", p).group(1))
        for p in glob.glob("ebg_p*-genos.txt")
    )
    if not groups:
        sys.exit("error: no EBG genotype files (ebg_p*-genos.txt) found")

    rows = []  # (snp, ind, ploidy, ebg_alt, updog_alt, status)
    for g in groups:
        samples = load_samples(f"p{g}.samples.txt")
        snps = load_loci(f"p{g}.loci.txt")
        ebg = load_ebg(f"ebg_p{g}-genos.txt", samples, snps)
        updog = load_updog(f"updog_p{g}.updog.tsv", g)
        for key in ebg:
            e = ebg[key]
            u = updog.get(key)
            if e is None or u is None:
                status = "missing"
            elif e == u:
                status = "agree"
            else:
                status = "disagree"
            rows.append((key[0], key[1], g, e, u, status))

    rows.sort(key=lambda r: (r[0], r[1]))
    with open("genotype_concordance.tsv", "w") as fh:
        fh.write("snp\tind\tploidy\tebg_alt_dosage\tupdog_alt_dosage\tstatus\n")
        for snp, ind, g, e, u, status in rows:
            fh.write(
                f"{snp}\t{ind}\t{g}\t"
                f"{'NA' if e is None else e}\t{'NA' if u is None else u}\t{status}\n"
            )

    # Summaries
    def rate(agree, total):
        return f"{agree / total:.4f}" if total else "NA"

    by_ploidy = defaultdict(lambda: [0, 0])       # ploidy -> [agree, compared]
    by_dosage = defaultdict(lambda: [0, 0])       # ebg_alt -> [agree, compared]
    discord_by_snp = defaultdict(int)
    total_agree = total_compared = total_missing = 0
    for snp, ind, g, e, u, status in rows:
        if status == "missing":
            total_missing += 1
            continue
        total_compared += 1
        agree = 1 if status == "agree" else 0
        total_agree += agree
        by_ploidy[g][0] += agree
        by_ploidy[g][1] += 1
        by_dosage[e][0] += agree
        by_dosage[e][1] += 1
        if status == "disagree":
            discord_by_snp[snp] += 1

    with open("genotype_concordance_summary.tsv", "w") as fh:
        fh.write("metric\tvalue\n")
        fh.write(f"calls_compared\t{total_compared}\n")
        fh.write(f"calls_missing\t{total_missing}\n")
        fh.write(f"calls_agree\t{total_agree}\n")
        fh.write(f"calls_disagree\t{total_compared - total_agree}\n")
        fh.write(f"overall_concordance\t{rate(total_agree, total_compared)}\n")
        for g in sorted(by_ploidy):
            a, t = by_ploidy[g]
            fh.write(f"concordance_ploidy_{g}\t{rate(a, t)}\n")
        for d in sorted(by_dosage):
            a, t = by_dosage[d]
            fh.write(f"concordance_ebg_alt_dosage_{d}\t{rate(a, t)}\n")
        fh.write("\n# most discordant SNPs (snp\tdisagreeing_calls)\n")
        for snp, n in sorted(discord_by_snp.items(), key=lambda kv: (-kv[1], kv[0]))[:20]:
            fh.write(f"{snp}\t{n}\n")

    sys.stderr.write(
        f"[compare_genotypes] {total_compared} calls compared, "
        f"{total_agree} agree, {total_compared - total_agree} disagree, "
        f"{total_missing} missing\n"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
