#!/usr/bin/env python3
"""
Normalize FASTQ headers for gzipped read files.

Typical use case:
- Input R1 header: @SRR17610144.1 1/1
- Output R1 header: @SRR17610144.1/1

- Input R2 header: @SRR17610144.1 2/2
- Output R2 header: @SRR17610144.1/2

For paired-end input, the script validates paired-end consistency by ensuring
paired records share an identical base read ID (ignoring trailing /1 or /2).
"""

from __future__ import annotations

import argparse
import gzip
import os
import re
import sys
from dataclasses import dataclass
from typing import Iterator, Optional, TextIO, Tuple


HEADER_SUFFIX_RE = re.compile(r"^(?P<base>.+?)(?:/(?P<mate>[12]))?$")
READ_HINT_RE = re.compile(r"^(?P<mate>[12])(?:$|[/:].*)")


@dataclass
class Stats:
    records: int = 0
    id_mismatches: int = 0
    r1_hint_conflicts: int = 0
    r2_hint_conflicts: int = 0
    rewritten_r1: int = 0
    rewritten_r2: int = 0


class FastqFormatError(RuntimeError):
    pass


def derive_output_path(input_path: str, mate: int) -> str:
    if input_path.endswith(".fastq.gz"):
        stem = input_path[: -len(".fastq.gz")]
    elif input_path.endswith(".fq.gz"):
        stem = input_path[: -len(".fq.gz")]
    elif input_path.endswith(".gz"):
        stem = input_path[: -len(".gz")]
    else:
        stem = input_path
    return f"{stem}.normalized_R{mate}.fastq.gz"


def parse_header(header_line: str) -> Tuple[str, Optional[str], str]:
    """
    Return (base_id, read_hint, first_token_without_at).

    - base_id: first token with optional trailing /1 or /2 removed
    - read_hint: parsed mate hint ('1'/'2') from either first or second token
    """
    line = header_line.rstrip("\n")
    if not line.startswith("@"):
        raise FastqFormatError(f"Header does not start with '@': {line}")

    fields = line[1:].split()
    if not fields:
        raise FastqFormatError("Empty FASTQ header line")

    first_token = fields[0]
    m = HEADER_SUFFIX_RE.match(first_token)
    if not m:
        raise FastqFormatError(f"Unable to parse FASTQ header token: {first_token}")

    base_id = m.group("base")
    hint = m.group("mate")

    if len(fields) > 1:
        m2 = READ_HINT_RE.match(fields[1])
        if m2 and hint is None:
            hint = m2.group("mate")

    return base_id, hint, first_token


def read_fastq_records(handle: TextIO, source_name: str) -> Iterator[Tuple[str, str, str, str, int]]:
    line_no = 0
    while True:
        h = handle.readline()
        if not h:
            return
        s = handle.readline()
        p = handle.readline()
        q = handle.readline()
        line_no += 4

        if not (s and p and q):
            raise FastqFormatError(
                f"Incomplete FASTQ record near line {line_no} in {source_name}"
            )
        if not p.startswith("+"):
            raise FastqFormatError(
                f"FASTQ third line does not start with '+' near line {line_no - 1} in {source_name}"
            )
        if len(s.rstrip("\n")) != len(q.rstrip("\n")):
            raise FastqFormatError(
                f"Sequence/quality length mismatch near line {line_no} in {source_name}"
            )

        yield h, s, p, q, line_no - 3


def normalize_pair(
    in_r1: str,
    in_r2: str,
    out_r1: str,
    out_r2: str,
    allow_id_mismatch: bool,
) -> Stats:
    stats = Stats()

    with gzip.open(in_r1, "rt") as r1_in, gzip.open(in_r2, "rt") as r2_in, gzip.open(
        out_r1, "wt"
    ) as r1_out, gzip.open(out_r2, "wt") as r2_out:
        it1 = read_fastq_records(r1_in, in_r1)
        it2 = read_fastq_records(r2_in, in_r2)

        while True:
            rec1 = next(it1, None)
            rec2 = next(it2, None)

            if rec1 is None and rec2 is None:
                break
            if rec1 is None or rec2 is None:
                raise FastqFormatError(
                    "R1 and R2 contain different numbers of FASTQ records"
                )

            h1, s1, p1, q1, line1 = rec1
            h2, s2, p2, q2, line2 = rec2

            base1, hint1, tok1 = parse_header(h1)
            base2, hint2, tok2 = parse_header(h2)

            if base1 != base2:
                stats.id_mismatches += 1
                msg = (
                    f"Pair ID mismatch at record {stats.records + 1}: "
                    f"R1={base1} (line {line1}) vs R2={base2} (line {line2})"
                )
                if not allow_id_mismatch:
                    raise FastqFormatError(msg)
                print(f"WARNING: {msg}", file=sys.stderr)

            if hint1 not in (None, "1"):
                stats.r1_hint_conflicts += 1
            if hint2 not in (None, "2"):
                stats.r2_hint_conflicts += 1

            new_h1 = f"@{base1}/1\n"
            new_h2 = f"@{base2}/2\n"

            if new_h1.rstrip("\n") != h1.rstrip("\n"):
                stats.rewritten_r1 += 1
            if new_h2.rstrip("\n") != h2.rstrip("\n"):
                stats.rewritten_r2 += 1

            r1_out.write(new_h1)
            r1_out.write(s1)
            r1_out.write(p1)
            r1_out.write(q1)

            r2_out.write(new_h2)
            r2_out.write(s2)
            r2_out.write(p2)
            r2_out.write(q2)

            stats.records += 1

    return stats


def normalize_single(in_r1: str, out_r1: str) -> Stats:
    stats = Stats()

    with gzip.open(in_r1, "rt") as r1_in, gzip.open(out_r1, "wt") as r1_out:
        for h1, s1, p1, q1, _line1 in read_fastq_records(r1_in, in_r1):
            base1, hint1, _tok1 = parse_header(h1)

            if hint1 not in (None, "1"):
                stats.r1_hint_conflicts += 1

            new_h1 = f"@{base1}/1\n"
            if new_h1.rstrip("\n") != h1.rstrip("\n"):
                stats.rewritten_r1 += 1

            r1_out.write(new_h1)
            r1_out.write(s1)
            r1_out.write(p1)
            r1_out.write(q1)
            stats.records += 1

    return stats


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Check FASTQ headers and rewrite them to @<base_id>/1 (single-end) or "
            "@<base_id>/1 and @<base_id>/2 (paired-end) for gzipped input FASTQs."
        )
    )
    parser.add_argument("--r1", required=True, help="Input R1 FASTQ (.gz)")
    parser.add_argument("--r2", default=None, help="Input R2 FASTQ (.gz)")
    parser.add_argument("--out-r1", default=None, help="Output R1 FASTQ (.gz)")
    parser.add_argument("--out-r2", default=None, help="Output R2 FASTQ (.gz)")
    parser.add_argument(
        "--allow-id-mismatch",
        action="store_true",
        help="Continue on R1/R2 ID mismatches (warn instead of failing)",
    )
    return parser


def validate_paths(args: argparse.Namespace) -> None:
    if not os.path.exists(args.r1):
        raise FileNotFoundError(f"R1 not found: {args.r1}")
    if args.r2 and not os.path.exists(args.r2):
        raise FileNotFoundError(f"R2 not found: {args.r2}")

    if args.out_r1 is None:
        args.out_r1 = derive_output_path(args.r1, 1)
    if args.r2 and args.out_r2 is None:
        args.out_r2 = derive_output_path(args.r2, 2)
    if not args.r2 and args.out_r2:
        raise ValueError("--out-r2 is only valid when --r2 is provided")

    output_paths = [args.out_r1]
    if args.out_r2:
        output_paths.append(args.out_r2)

    for p in output_paths:
        out_dir = os.path.dirname(os.path.abspath(p))
        if out_dir and not os.path.isdir(out_dir):
            raise FileNotFoundError(f"Output directory not found: {out_dir}")


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    try:
        validate_paths(args)
        if args.r2:
            stats = normalize_pair(
                in_r1=args.r1,
                in_r2=args.r2,
                out_r1=args.out_r1,
                out_r2=args.out_r2,
                allow_id_mismatch=args.allow_id_mismatch,
            )
        else:
            stats = normalize_single(
                in_r1=args.r1,
                out_r1=args.out_r1,
            )
    except Exception as exc:  # noqa: BLE001
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print("FASTQ header normalization complete", file=sys.stderr)
    print(f"Records processed: {stats.records}", file=sys.stderr)
    print(f"R1 headers rewritten: {stats.rewritten_r1}", file=sys.stderr)
    print(f"R2 headers rewritten: {stats.rewritten_r2}", file=sys.stderr)
    print(f"R1 hint conflicts: {stats.r1_hint_conflicts}", file=sys.stderr)
    print(f"R2 hint conflicts: {stats.r2_hint_conflicts}", file=sys.stderr)
    print(f"R1/R2 ID mismatches: {stats.id_mismatches}", file=sys.stderr)
    print(f"Output R1: {args.out_r1}", file=sys.stderr)
    if args.out_r2:
        print(f"Output R2: {args.out_r2}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
