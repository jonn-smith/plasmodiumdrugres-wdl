#!/usr/bin/env python3
"""Self-check for the sample-table plumbing in resmarkers_to_terra.py.

Run with: python3 python/test_resmarkers_to_terra.py
No pytest, no fixtures, no ops_utils needed (the Terra import is deferred).
"""

import csv
import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import resmarkers_to_terra as r

CALLS = (
    "specimen_name\tgene\taa_position\taa\n"
    "SPEC-1\tdhfr-ts\t164\tL\n"
    "SPEC-1\tdhfr-ts\t164\tI\n"   # mixed call for the same gene/position
    "SPEC-2\tmdr1\t86\tY\n"
)


# The real row ids from the mad4hatter_sample table of
# broad-firecloud-dsde-methods/JTS_COPY_20261024_..._Mad4Hatter_CLONE, where upserting
# the bare specimen_name created five half-empty duplicate rows.
REAL_TERRA_IDS = {
    "SRR26819135_S1_L001", "SRR26819139_S2_L001", "SRR26819141_S3_L001",
    "SRR26819151_S4_L001", "SRR26819200_S5_L001",
}


def test_resolve_terra_ids():
    # the regression: a bare specimen_name must resolve onto the suffixed row id
    resolved = r.resolve_terra_ids(
        {"SRR26819135", "SRR26819139", "SRR26819141", "SRR26819151", "SRR26819200"},
        REAL_TERRA_IDS)
    assert resolved == {
        "SRR26819135": "SRR26819135_S1_L001",
        "SRR26819139": "SRR26819139_S2_L001",
        "SRR26819141": "SRR26819141_S3_L001",
        "SRR26819151": "SRR26819151_S4_L001",
        "SRR26819200": "SRR26819200_S5_L001",
    }, resolved

    # an id that already matches exactly is used as-is
    assert r.resolve_terra_ids({"SRR26819135_S1_L001"}, REAL_TERRA_IDS) == {
        "SRR26819135_S1_L001": "SRR26819135_S1_L001"}

    # lane-less suffix also matches
    assert r.resolve_terra_ids({"S"}, {"S_S7"}) == {"S": "S_S7"}

    # ambiguous and unmatched names are skipped, never guessed
    assert r.resolve_terra_ids({"X"}, {"X_S1", "X_S2"}) == {}
    assert r.resolve_terra_ids({"NOPE"}, REAL_TERRA_IDS) == {}

    # a near-miss must NOT match: suffix has to be _S<digits>
    assert r.resolve_terra_ids({"SRR26819135"}, {"SRR26819135_extra"}) == {}


def main():
    test_resolve_terra_ids()

    # Terra demands the row id column be named <table>_id; ops_utils raises otherwise.
    assert r.id_column_for("sample") == "sample_id"
    assert r.id_column_for("specimen") == "specimen_id"

    with tempfile.TemporaryDirectory() as d:
        calls_tsv = os.path.join(d, "collapsed_amino_acid_calls.tsv")
        with open(calls_tsv, "w") as fh:
            fh.write(CALLS)

        for table in ("sample", "specimen"):
            id_column = r.id_column_for(table)
            rows, columns = r.pivot(calls_tsv, id_column)

            # the row id lands under the detected column, not a hardcoded "sample_id"
            assert sorted(row[id_column] for row in rows) == ["SPEC-1", "SPEC-2"], rows
            assert columns == ["resmarker_dhfr-ts_164", "resmarker_mdr1_86"], columns
            spec1 = next(row for row in rows if row[id_column] == "SPEC-1")
            assert spec1["resmarker_dhfr-ts_164"] == "I/L", spec1

            out_tsv = os.path.join(d, f"{table}.tsv")
            r.write_tsv(out_tsv, rows, columns, id_column)
            with open(out_tsv) as fh:
                header = next(csv.reader(fh, delimiter="\t"))
            assert header[0] == id_column, header

    print("OK")


if __name__ == "__main__":
    main()
