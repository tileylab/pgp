#!/usr/bin/env python3
"""
Generate a minimal synthetic test dataset for the PGP pipeline.

Creates:
  - reference.fa        : 3 contigs (~500bp each)
  - intervals.list      : contig names
  - samplesheet.csv     : 2 paired-end samples
  - sample1_R{1,2}.fastq.gz, sample2_R{1,2}.fastq.gz

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

# SNP positions per sample (0-indexed within contig)
# sample1 has SNPs on chr1; sample2 has SNPs on chr1 and chr2
SAMPLE_SNPS = {
    "sample1": {"chr1": {100: "T", 250: "G"}},
    "sample2": {"chr1": {100: "T", 350: "A"}, "chr2": {200: "C"}},
}


def random_seq(length):
    return "".join(random.choice("ACGT") for _ in range(length))


def reverse_complement(seq):
    comp = str.maketrans("ACGT", "TGCA")
    return seq.translate(comp)[::-1]


def mutate_seq(seq, contig, snps_dict):
    """Apply SNPs for a given contig."""
    seq_list = list(seq)
    if contig in snps_dict:
        for pos, alt in snps_dict[contig].items():
            if pos < len(seq_list):
                seq_list[pos] = alt
    return "".join(seq_list)


def make_quality(length):
    """Generate a plausible quality string (Phred+33, mostly high quality)."""
    return "".join(chr(random.randint(53, 73)) for _ in range(length))


def generate_reads(ref_seqs, sample_snps, n_reads):
    """Generate paired-end reads from reference with optional SNPs."""
    r1_records = []
    r2_records = []
    contigs = list(ref_seqs.keys())
    contig_lengths = [len(ref_seqs[c]) for c in contigs]
    total_len = sum(contig_lengths)

    for i in range(n_reads):
        # Pick contig weighted by length
        r = random.random() * total_len
        cumulative = 0
        chosen_contig = contigs[0]
        for c, cl in zip(contigs, contig_lengths):
            cumulative += cl
            if r < cumulative:
                chosen_contig = c
                break

        seq = ref_seqs[chosen_contig]
        # Apply sample-specific SNPs
        seq = mutate_seq(seq, chosen_contig, sample_snps)

        max_start = len(seq) - FRAG_LEN
        if max_start < 0:
            max_start = 0
        start = random.randint(0, max_start)
        fragment = seq[start : start + FRAG_LEN]

        # R1 is forward from fragment start, R2 is reverse complement from fragment end
        r1_seq = fragment[:READ_LEN]
        r2_seq = reverse_complement(fragment[-READ_LEN:])

        read_name = f"@read_{i}/1"
        r1_records.append(f"{read_name}\n{r1_seq}\n+\n{make_quality(READ_LEN)}")

        read_name = f"@read_{i}/2"
        r2_records.append(f"{read_name}\n{r2_seq}\n+\n{make_quality(READ_LEN)}")

    return r1_records, r2_records


def write_fastq_gz(records, filepath):
    with gzip.open(filepath, "wt") as f:
        f.write("\n".join(records) + "\n")


def main():
    # Generate reference
    ref_seqs = {}
    ref_lines = []
    interval_lines = []

    for contig, length in CONTIGS.items():
        seq = random_seq(length)
        ref_seqs[contig] = seq
        ref_lines.append(f">{contig}")
        # Write sequence in 80-char lines
        for j in range(0, len(seq), 80):
            ref_lines.append(seq[j : j + 80])
        interval_lines.append(contig)

    ref_path = os.path.join(SCRIPT_DIR, "reference.fa")
    with open(ref_path, "w") as f:
        f.write("\n".join(ref_lines) + "\n")
    print(f"Wrote {ref_path}")

    intervals_path = os.path.join(SCRIPT_DIR, "intervals.list")
    with open(intervals_path, "w") as f:
        f.write("\n".join(interval_lines) + "\n")
    print(f"Wrote {intervals_path}")

    # Generate reads for each sample
    for sample, snps in SAMPLE_SNPS.items():
        r1, r2 = generate_reads(ref_seqs, snps, N_READS_PER_SAMPLE)
        r1_path = os.path.join(SCRIPT_DIR, f"{sample}_R1.fastq.gz")
        r2_path = os.path.join(SCRIPT_DIR, f"{sample}_R2.fastq.gz")
        write_fastq_gz(r1, r1_path)
        write_fastq_gz(r2, r2_path)
        print(f"Wrote {r1_path}, {r2_path}")

    # Generate samplesheet (relative paths — nf-schema resolves them relative to samplesheet location)
    samplesheet_path = os.path.join(SCRIPT_DIR, "samplesheet.csv")
    # sample2 is flagged as an outgroup (1) to exercise the ingroup-only products
    outgroup_flags = {"sample1": 0, "sample2": 1}
    with open(samplesheet_path, "w") as f:
        f.write("sample,fastq_1,fastq_2,species,ploidy,outgroup\n")
        for sample in SAMPLE_SNPS:
            og = outgroup_flags.get(sample, 0)
            f.write(f"{sample},{sample}_R1.fastq.gz,{sample}_R2.fastq.gz,test_species,2,{og}\n")
    print(f"Wrote {samplesheet_path}")


if __name__ == "__main__":
    main()
