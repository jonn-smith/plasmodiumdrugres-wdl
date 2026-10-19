#!/usr/bin/env python3
"""Reformat a MAD4HatTeR allele table (e.g. allele_data_filtered.txt)
into the specimen/target/reads/seq TSV layout.

Input  (tab-separated): sample_name, target_name, reads, <seq column>, ...
Output (tab-separated): specimen_name, target_name, reads, seq

Columns are renamed as below, reordered, and all other columns are dropped.
Values are copied through unchanged.

The seq column defaults to the unmasked "asv" column. Masked sequences
(asv_masked, with lowercase 'n') break translate_loci_of_interest.R, whose
alignment only scores A/C/G/T.

Only the Python standard library is used.
"""

import argparse
import csv
import sys


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Reformat an allele table into specimen_name/target_name/reads/seq TSV."
    )
    parser.add_argument("--input", required=True,
                        help="Tab-separated allele table (sample_name, target_name, reads, asv, ...)")
    parser.add_argument("--output", required=True, help="Output TSV path")
    parser.add_argument("--seq-column", default="asv",
                        help="Input column to write as seq (default: asv, the unmasked sequence)")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    # output column -> input column, in output order
    column_map = (
        ("specimen_name", "sample_name"),
        ("target_name", "target_name"),
        ("reads", "reads"),
        ("seq", args.seq_column),
    )
    with open(args.input, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = [src for _, src in column_map if src not in (reader.fieldnames or [])]
        if missing:
            sys.exit(f"ERROR: {args.input} is missing required column(s): {', '.join(missing)}")

        with open(args.output, "w", newline="") as out:
            writer = csv.writer(out, delimiter="\t", lineterminator="\n")
            writer.writerow([dst for dst, _ in column_map])
            n = 0
            for row in reader:
                writer.writerow([row[src] for _, src in column_map])
                n += 1
    print(f"Wrote {n} rows to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
