#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import re
import csv
import sys

FILE = Path("Scientific_Data_manuscript_working.txt")
OUTDIR = Path("18B_reference_citation_audit")


def fail(msg):
    print("ERROR:", msg, file=sys.stderr)
    sys.exit(1)


if not FILE.is_file():
    fail(f"Manuscript not found: {FILE}")

text = FILE.read_text(encoding="utf-8")

if "\nREFERENCES\n" not in text:
    fail("REFERENCES heading not found.")

body, refs_and_after = text.split("\nREFERENCES\n", 1)

if "\nFIGURE LEGENDS\n" not in refs_and_after:
    fail("FIGURE LEGENDS heading not found.")

refs, legends = refs_and_after.split(
    "\nFIGURE LEGENDS\n",
    1
)


# ============================================================
# 1. Normalize Markdown DOI links in references only
#
# [https://doi.org/xxx](https://doi.org/xxx)
# ->
# https://doi.org/xxx
# ============================================================

refs_clean = re.sub(
    r"\[(https://doi\.org/[^\]]+)\]"
    r"\(\1\)",
    r"\1",
    refs,
    flags=re.I
)


# ============================================================
# 2. Extract reference list
# ============================================================

reference_entries = {}

for line in refs_clean.splitlines():

    line = line.strip()

    m = re.match(
        r"^(\d+)\.\s+(.+)$",
        line
    )

    if m:
        n = int(m.group(1))
        reference_entries[n] = m.group(2)


if not reference_entries:
    fail("No numbered references detected.")

max_ref = max(reference_entries)

expected_refs = list(
    range(1, max_ref + 1)
)

observed_refs = sorted(
    reference_entries
)

references_continuous = (
    observed_refs == expected_refs
)


# ============================================================
# 3. Expand citation groups
#
# Supports:
# [1]
# [1,2]
# [1, 2, 5]
# [1–3]
# [1-3]
# [1–3,5]
# ============================================================

def expand_citation_group(content):

    content = (
        content
        .replace("–", "-")
        .replace("—", "-")
    )

    nums = []

    for item in re.split(
        r"[,;]",
        content
    ):

        item = item.strip()

        if not item:
            continue

        m_range = re.fullmatch(
            r"(\d+)\s*-\s*(\d+)",
            item
        )

        if m_range:

            start = int(
                m_range.group(1)
            )

            end = int(
                m_range.group(2)
            )

            if end < start:
                raise ValueError(
                    f"Reverse citation range: {item}"
                )

            nums.extend(
                range(start, end + 1)
            )

            continue

        if re.fullmatch(
            r"\d+",
            item
        ):
            nums.append(
                int(item)
            )
            continue

        raise ValueError(
            f"Cannot parse citation item: {item}"
        )

    return nums


citation_pattern = re.compile(
    r"\["
    r"(\s*\d+"
    r"(?:\s*[-–—]\s*\d+)?"
    r"(?:\s*[,;]\s*\d+"
    r"(?:\s*[-–—]\s*\d+)?)*"
    r"\s*)"
    r"\]"
)


citation_occurrences = []

for m in citation_pattern.finditer(body):

    raw = m.group(1)

    try:
        nums = expand_citation_group(
            raw
        )
    except ValueError as e:
        fail(str(e))

    citation_occurrences.append({
        "position": m.start(),
        "raw": "[" + raw + "]",
        "numbers": nums,
    })


all_cited = []

for x in citation_occurrences:
    all_cited.extend(
        x["numbers"]
    )

cited_set = set(
    all_cited
)


# ============================================================
# 4. Undefined / uncited references
# ============================================================

undefined = sorted(
    x
    for x in cited_set
    if x not in reference_entries
)

uncited = sorted(
    x
    for x in reference_entries
    if x not in cited_set
)


# ============================================================
# 5. First-citation order
# ============================================================

first_seen = []

seen = set()

for occurrence in citation_occurrences:

    for n in occurrence["numbers"]:

        if n not in seen:
            seen.add(n)
            first_seen.append(n)


expected_first_order = [
    x
    for x in expected_refs
    if x in cited_set
]

first_order_ok = (
    first_seen == expected_first_order
)


# ============================================================
# 6. DOI checks
# ============================================================

doi_pattern = re.compile(
    r"https://doi\.org/"
    r"(10\.\d{4,9}/[-._;()/:A-Za-z0-9]+)",
    flags=re.I
)

doi_rows = []

doi_to_refs = {}

for n, entry in reference_entries.items():

    dois = doi_pattern.findall(
        entry
    )

    doi_rows.append({
        "reference": n,
        "doi_count": len(dois),
        "dois": ";".join(dois),
    })

    for doi in dois:

        key = doi.lower().rstrip(".")

        doi_to_refs.setdefault(
            key,
            []
        ).append(n)


duplicate_dois = {
    doi: refs_
    for doi, refs_ in doi_to_refs.items()
    if len(refs_) > 1
}


# ============================================================
# 7. Year checks
# ============================================================

missing_year = []

for n, entry in reference_entries.items():

    if not re.search(
        r"\((?:19|20)\d{2}\)",
        entry
    ):
        missing_year.append(n)


# ============================================================
# 8. Write cleaned manuscript
# ============================================================

new_text = (
    body.rstrip()
    + "\n\nREFERENCES\n\n"
    + refs_clean.strip()
    + "\n\nFIGURE LEGENDS\n\n"
    + legends.lstrip()
)

FILE.write_text(
    new_text,
    encoding="utf-8"
)


# ============================================================
# 9. Reports
# ============================================================

OUTDIR.mkdir(
    parents=True,
    exist_ok=True
)

report = []


def add(check, status, observed, expected):

    report.append({
        "check": check,
        "status": status,
        "observed": str(observed),
        "expected": str(expected),
    })


add(
    "reference_numbering_continuous",
    "PASS" if references_continuous else "FAIL",
    observed_refs,
    expected_refs
)

add(
    "undefined_citations",
    "PASS" if not undefined else "FAIL",
    undefined,
    []
)

add(
    "uncited_references",
    "PASS" if not uncited else "FAIL",
    uncited,
    []
)

add(
    "first_citation_order",
    "PASS" if first_order_ok else "FAIL",
    first_seen,
    expected_first_order
)

add(
    "duplicate_DOIs",
    "PASS" if not duplicate_dois else "FAIL",
    duplicate_dois,
    {}
)

add(
    "reference_years",
    "PASS" if not missing_year else "WARN",
    missing_year,
    []
)


report_file = (
    OUTDIR
    / "18B_reference_citation_report.tsv"
)

with report_file.open(
    "w",
    encoding="utf-8",
    newline=""
) as f:

    writer = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=[
            "check",
            "status",
            "observed",
            "expected",
        ],
        lineterminator="\n"
    )

    writer.writeheader()
    writer.writerows(report)


occ_file = (
    OUTDIR
    / "18B_citation_occurrences.tsv"
)

with occ_file.open(
    "w",
    encoding="utf-8",
    newline=""
) as f:

    writer = csv.writer(
        f,
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writerow([
        "order",
        "raw_citation",
        "expanded_references",
    ])

    for i, x in enumerate(
        citation_occurrences,
        1
    ):

        writer.writerow([
            i,
            x["raw"],
            ",".join(
                map(str, x["numbers"])
            ),
        ])


doi_file = (
    OUTDIR
    / "18B_reference_DOI_audit.tsv"
)

with doi_file.open(
    "w",
    encoding="utf-8",
    newline=""
) as f:

    writer = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=[
            "reference",
            "doi_count",
            "dois",
        ],
        lineterminator="\n"
    )

    writer.writeheader()
    writer.writerows(
        doi_rows
    )


P = sum(
    x["status"] == "PASS"
    for x in report
)

W = sum(
    x["status"] == "WARN"
    for x in report
)

F = sum(
    x["status"] == "FAIL"
    for x in report
)

overall = (
    "FAIL"
    if F
    else (
        "WARN"
        if W
        else "PASS"
    )
)


log_file = (
    OUTDIR
    / "18B_reference_citation_audit.log"
)

with log_file.open(
    "w",
    encoding="utf-8"
) as f:

    f.write(
        "=" * 78 + "\n"
    )

    f.write(
        "FINAL REFERENCE & CITATION AUDIT\n"
    )

    f.write(
        "=" * 78 + "\n"
    )

    f.write(
        f"OVERALL: {overall}\n"
    )

    f.write(
        f"PASS={P} WARN={W} FAIL={F}\n"
    )

    f.write(
        f"Reference count: {len(reference_entries)}\n"
    )

    f.write(
        f"Cited references: {sorted(cited_set)}\n"
    )

    f.write(
        f"First-use order: {first_seen}\n"
    )

    for x in report:

        if x["status"] != "PASS":

            f.write(
                f"[{x['status']}] "
                f"{x['check']}: "
                f"observed={x['observed']} "
                f"expected={x['expected']}\n"
            )


print("=" * 78)
print("FINAL REFERENCE & CITATION AUDIT")
print("=" * 78)

print(
    f"OVERALL: {overall}   "
    f"PASS={P} WARN={W} FAIL={F}"
)

print(
    "Reference numbers :",
    observed_refs
)

print(
    "Cited references  :",
    sorted(cited_set)
)

print(
    "First-use order   :",
    first_seen
)

print(
    "Undefined         :",
    undefined
)

print(
    "Uncited           :",
    uncited
)

print(
    "Duplicate DOIs    :",
    duplicate_dois
)

print(
    "Missing year      :",
    missing_year
)

print("-" * 78)
print("Report :", report_file)
print("DOI    :", doi_file)
print("Log    :", log_file)
print("=" * 78)

sys.exit(
    1 if F else 0
)
