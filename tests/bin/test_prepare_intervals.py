#!/usr/bin/env python3
"""Tests for bin/prepare_intervals.py"""

import os
import subprocess
import sys
import tempfile

SCRIPT = os.path.join(
    os.path.dirname(__file__), "..", "..", "bin", "prepare_intervals.py"
)


def run_prepare_intervals(fasta_content):
    """Run prepare_intervals.py on a temp FASTA and return output files."""
    with tempfile.TemporaryDirectory() as tmpdir:
        fasta_path = os.path.join(tmpdir, "input.fa")
        with open(fasta_path, "w") as f:
            f.write(fasta_content)

        subprocess.run(
            [sys.executable, SCRIPT, fasta_path],
            check=True,
            cwd=tmpdir,
        )

        results = {}
        for fname in ["intervals.list", "interval_map.txt", "ref.fa"]:
            path = os.path.join(tmpdir, fname)
            if os.path.exists(path):
                with open(path) as f:
                    results[fname] = f.read()
        return results


def test_basic_two_contigs():
    """Two contigs should produce numeric renaming and correct intervals."""
    fasta = ">contigA\nACGTACGT\n>contigB\nTTTTAAAA\n"
    results = run_prepare_intervals(fasta)

    # intervals.list should have two numeric entries
    intervals = results["intervals.list"].strip().split("\n")
    assert intervals == ["1", "2"]

    # interval_map should map numbers to original names
    mapping = results["interval_map.txt"].strip().split("\n")
    assert mapping == ["1\tcontigA", "2\tcontigB"]

    # ref.fa should have numeric headers with same sequences
    lines = results["ref.fa"].strip().split("\n")
    assert lines[0] == ">1"
    assert lines[1] == "ACGTACGT"
    assert lines[2] == ">2"
    assert lines[3] == "TTTTAAAA"


def test_single_contig():
    """Single contig input."""
    fasta = ">chr1\nGGGGCCCC\n"
    results = run_prepare_intervals(fasta)

    intervals = results["intervals.list"].strip().split("\n")
    assert intervals == ["1"]

    lines = results["ref.fa"].strip().split("\n")
    assert lines[0] == ">1"
    assert lines[1] == "GGGGCCCC"


def test_multiline_sequence():
    """Sequence split across multiple lines should be preserved."""
    fasta = ">seqX\nAAAA\nCCCC\nGGGG\n"
    results = run_prepare_intervals(fasta)

    ref_lines = results["ref.fa"].strip().split("\n")
    assert ref_lines[0] == ">1"
    # Each line is written separately (script writes line by line)
    assert ref_lines[1] == "AAAA"
    assert ref_lines[2] == "CCCC"
    assert ref_lines[3] == "GGGG"


def test_complex_headers():
    """Headers with spaces and descriptions should be fully captured."""
    fasta = ">chr1 some description here\nACGT\n>chr2 another desc\nTGCA\n"
    results = run_prepare_intervals(fasta)

    mapping = results["interval_map.txt"].strip().split("\n")
    assert mapping[0] == "1\tchr1 some description here"
    assert mapping[1] == "2\tchr2 another desc"


def test_all_output_files_created():
    """All three output files should be created."""
    fasta = ">a\nACGT\n"
    results = run_prepare_intervals(fasta)
    assert "intervals.list" in results
    assert "interval_map.txt" in results
    assert "ref.fa" in results
