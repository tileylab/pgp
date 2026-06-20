#!/usr/bin/env python3
"""
Generate a minimal synthetic test dataset for the PGP pipeline.

Creates:
  - reference.fa        : 3 contigs
  - intervals.list      : contig names
  - samplesheet.csv     : 6 paired-end samples (sample1-5 ingroup ploidy 4,
                          sample6 outgroup ploidy 2)
  - sample{1..6}_R{1,2}.fastq.gz

Reads carry real **dosage signal** at a shared set of biallelic SNP sites: for
each sample each SNP is assigned an ALT-allele dosage (0..ploidy), and every read
covering that site independently shows the ALT base with probability
dosage/ploidy. The resulting VCF `AD` field therefore reflects the dosage, which
is what EBG and updog genotype from. tests/default.nf.test only uses sample1 and
sample2 (fast); tests/genotyping.nf.test uses the full cohort.

Usage:
    python generate_test_data.py
"""

import gzip
import os
import random

random.seed(42)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
READ_LEN = 150
FRAG_LEN = 300
N_READS_PER_SAMPLE = 500
CONTIGS = {
    "chr1": 600,
    "chr2": 500,
    "chr3": 400,
}

# Samples: sample1-5 ingroup tetraploid, sample6 outgroup diploid.
SAMPLES = {
    "sample1": {"ploidy": 4, "outgroup": 0},
    "sample2": {"ploidy": 4, "outgroup": 0},
    "sample3": {"ploidy": 4, "outgroup": 0},
    "sample4": {"ploidy": 4, "outgroup": 0},
    "sample5": {"ploidy": 4, "outgroup": 0},
    "sample6": {"ploidy": 2, "outgroup": 1},
}

# Shared biallelic SNP sites (0-indexed positions per contig).
SNP_POSITIONS = {
    "chr1": [100, 250, 400],
    "chr2": [150, 300],
    "chr3": [200],
}


def random_seq(length):
    return "".join(random.choice("ACGT") for _ in range(length))


def reverse_complement(seq):
    comp = str.maketrans("ACGT", "TGCA")
    return seq.translate(comp)[::-1]


def make_quality(length):
    """Plausible high-quality Phred+33 string."""
    return "".join(chr(random.randint(53, 73)) for _ in range(length))


def generate_reads(ref_seqs, alt_bases, dosage, ploidy, n_reads):
    """Paired-end reads with per-read ALT draws at the SNP sites.

    alt_bases : {(contig, pos): alt_base}
    dosage    : {(contig, pos): alt-allele dosage 0..ploidy} for this sample
    """
    r1_records, r2_records = [], []
    contigs = list(ref_seqs.keys())
    contig_lengths = [len(ref_seqs[c]) for c in contigs]
    total_len = sum(contig_lengths)

    for i in range(n_reads):
        r = random.random() * total_len
        cumulative = 0
        chosen = contigs[0]
        for c, cl in zip(contigs, contig_lengths):
            cumulative += cl
            if r < cumulative:
                chosen = c
                break

        seq = ref_seqs[chosen]
        max_start = max(len(seq) - FRAG_LEN, 0)
        start = random.randint(0, max_start)
        fragment = list(seq[start : start + FRAG_LEN])

        # Independent ALT draw per covered SNP, at rate dosage/ploidy.
        for (c, pos), alt in alt_bases.items():
            if c == chosen and start <= pos < start + len(fragment):
                if random.random() < dosage[(c, pos)] / ploidy:
                    fragment[pos - start] = alt
        fragment = "".join(fragment)

        r1_seq = fragment[:READ_LEN]
        r2_seq = reverse_complement(fragment[-READ_LEN:])
        r1_records.append(f"@read_{i}/1\n{r1_seq}\n+\n{make_quality(READ_LEN)}")
        r2_records.append(f"@read_{i}/2\n{r2_seq}\n+\n{make_quality(READ_LEN)}")

    return r1_records, r2_records


def write_fastq_gz(records, filepath):
    with gzip.open(filepath, "wt") as f:
        f.write("\n".join(records) + "\n")


def main():
    # Reference
    ref_seqs, ref_lines, interval_lines = {}, [], []
    for contig, length in CONTIGS.items():
        seq = random_seq(length)
        ref_seqs[contig] = seq
        ref_lines.append(f">{contig}")
        for j in range(0, len(seq), 80):
            ref_lines.append(seq[j : j + 80])
        interval_lines.append(contig)

    with open(os.path.join(SCRIPT_DIR, "reference.fa"), "w") as f:
        f.write("\n".join(ref_lines) + "\n")
    with open(os.path.join(SCRIPT_DIR, "intervals.list"), "w") as f:
        f.write("\n".join(interval_lines) + "\n")

    # ALT base per SNP site (differs from the reference base there).
    alt_bases = {}
    for contig, positions in SNP_POSITIONS.items():
        for pos in positions:
            ref_base = ref_seqs[contig][pos]
            alt_bases[(contig, pos)] = random.choice([b for b in "ACGT" if b != ref_base])

    # Per-sample ALT dosage at each SNP, kept non-monomorphic across the ingroup.
    site_keys = list(alt_bases.keys())
    dosages = {s: {} for s in SAMPLES}
    for key in site_keys:
        while True:
            for s, meta in SAMPLES.items():
                dosages[s][key] = random.randint(0, meta["ploidy"])
            ingroup_vals = {dosages[s][key] for s, m in SAMPLES.items() if m["outgroup"] == 0}
            if len(ingroup_vals) > 1:  # avoid a monomorphic (uninformative) site
                break

    # Reads
    for sample, meta in SAMPLES.items():
        r1, r2 = generate_reads(ref_seqs, alt_bases, dosages[sample], meta["ploidy"], N_READS_PER_SAMPLE)
        write_fastq_gz(r1, os.path.join(SCRIPT_DIR, f"{sample}_R1.fastq.gz"))
        write_fastq_gz(r2, os.path.join(SCRIPT_DIR, f"{sample}_R2.fastq.gz"))
        print(f"Wrote {sample} reads (ploidy {meta['ploidy']}, outgroup {meta['outgroup']})")

    # Samplesheet
    with open(os.path.join(SCRIPT_DIR, "samplesheet.csv"), "w") as f:
        f.write("sample,fastq_1,fastq_2,species,ploidy,outgroup\n")
        for sample, meta in SAMPLES.items():
            f.write(
                f"{sample},{sample}_R1.fastq.gz,{sample}_R2.fastq.gz,"
                f"test_species,{meta['ploidy']},{meta['outgroup']}\n"
            )
    print("Wrote reference.fa, intervals.list, samplesheet.csv")


if __name__ == "__main__":
    main()
