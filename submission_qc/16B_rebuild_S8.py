#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import csv
import shutil
from datetime import datetime


src = Path(
    "submission_package_scientific_data/"
    "03_Supplementary_Tables/"
    "Supplementary_Table_S8_"
    "sample_exclusion_and_leave_one_out_summary.tsv"
)

if not src.exists():
    raise SystemExit(f"ERROR: file not found: {src}")


# ============================================================
# Backup the malformed version
# ============================================================

backup = src.with_suffix(
    src.suffix + ".pre16B.bak"
)

if not backup.exists():
    shutil.copy2(src, backup)
    print("Backup:", backup)
else:
    print("Backup already exists:", backup)


# ============================================================
# Canonical S8 structure
# ============================================================

header = [
    "Analysis_type",
    "Scenario_or_removed_sample",
    "N_samples",
    "PERMANOVA_R2",
    "PERMANOVA_P",
    "PERMDISP_P",
    "Cross_cohort_matrix_correlation",
    "Cross_cohort_matrix_MAE",
    "Minimum_genomewide_LFC_Spearman",
    "Review_type",
]


rows = [

    # --------------------------------------------------------
    # Primary / multi-sample sensitivity
    # --------------------------------------------------------

    [
        "Multi-sample exclusion",
        "Primary_all77",
        77,
        0.415491,
        "1e-04",
        0.7220,
        1.0,
        0.0,
        "NA",
        "Primary",
    ],

    [
        "Multi-sample exclusion",
        "Exclude_multi_metric",
        72,
        0.466386,
        "1e-04",
        0.0575,
        0.961368,
        0.054870,
        0.14051,
        "Expression review",
    ],

    [
        "Multi-sample exclusion",
        "Exclude_R2h4",
        76,
        0.416445,
        "1e-04",
        0.7371,
        0.999965,
        0.000433,
        0.00606,
        "Sequencing review",
    ],

    [
        "Multi-sample exclusion",
        "Exclude_union6",
        71,
        0.468294,
        "1e-04",
        0.0646,
        0.961648,
        0.054600,
        0.14054,
        "Expression + sequencing review",
    ],

    # --------------------------------------------------------
    # Single-sample leave-one-out
    # --------------------------------------------------------

    [
        "Single-sample leave-one-out",
        "RB_2",
        76,
        0.432401,
        "1e-04",
        0.3466,
        0.978515,
        0.06068,
        0.89220,
        "Expression review",
    ],

    [
        "Single-sample leave-one-out",
        "R0h_6",
        76,
        0.426524,
        "1e-04",
        0.5132,
        0.999020,
        0.00631,
        0.95450,
        "Expression review",
    ],

    [
        "Single-sample leave-one-out",
        "R6h_8",
        76,
        0.420773,
        "1e-04",
        0.7808,
        0.995241,
        0.00526,
        0.93773,
        "Expression review",
    ],

    [
        "Single-sample leave-one-out",
        "R10h_1",
        76,
        0.417850,
        "1e-04",
        0.6362,
        0.999230,
        0.00235,
        0.96378,
        "Expression review",
    ],

    [
        "Single-sample leave-one-out",
        "R18h_5",
        76,
        0.425804,
        "1e-04",
        0.5821,
        0.987333,
        0.00534,
        0.90111,
        "Expression review",
    ],

    [
        "Single-sample leave-one-out",
        "R2h_4",
        76,
        0.416445,
        "1e-04",
        0.7371,
        0.999965,
        0.000433,
        0.98733,
        "Sequencing review",
    ],
]


# ============================================================
# Internal checks before writing
# ============================================================

assert len(rows) == 10

keys = [
    (r[0], r[1])
    for r in rows
]

assert len(keys) == len(set(keys))

for r in rows:

    n = int(r[2])
    r2 = float(r[3])
    pp = float(r[4])
    pd = float(r[5])
    mc = float(r[6])
    mae = float(r[7])

    assert 1 <= n <= 77
    assert 0 <= r2 <= 1
    assert 0 <= pp <= 1
    assert 0 <= pd <= 1
    assert -1 <= mc <= 1
    assert mae >= 0


# ============================================================
# Write clean TSV
# ============================================================

with src.open(
    "w",
    encoding="utf-8",
    newline=""
) as f:

    w = csv.writer(
        f,
        delimiter="\t",
        lineterminator="\n"
    )

    w.writerow(header)
    w.writerows(rows)


# ============================================================
# Read-back QC
# ============================================================

with src.open(
    "r",
    encoding="utf-8",
    newline=""
) as f:

    reader = csv.DictReader(
        f,
        delimiter="\t"
    )

    check = list(reader)

assert len(check) == 10
assert reader.fieldnames == header

print("=" * 72)
print("S8 REBUILD COMPLETED")
print("=" * 72)
print("Output :", src)
print("Backup :", backup)
print("Rows   :", len(check))
print("Columns:", len(header))
print("STATUS : PASS")
print("=" * 72)
