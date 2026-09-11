#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import re
import sys
import zipfile
from pathlib import Path
from collections import defaultdict


# ============================================================
# Helpers
# ============================================================

def normalize(text):
    text = text.replace("\r\n", "\n").replace("\r", "\n")
    text = text.replace("\u00a0", " ")
    text = text.replace("–", "-")
    text = text.replace("—", "-")
    return text


def read_docx(path):

    parts = []

    with zipfile.ZipFile(path, "r") as z:

        for name in [
            "word/document.xml",
            "word/footnotes.xml",
            "word/endnotes.xml",
        ]:

            if name not in z.namelist():
                continue

            xml = z.read(name).decode(
                "utf-8",
                errors="replace"
            )

            xml = re.sub(
                r"</w:p>",
                "\n",
                xml
            )

            xml = re.sub(
                r"</w:tr>",
                "\n",
                xml
            )

            xml = re.sub(
                r"<w:tab[^>]*/>",
                "\t",
                xml
            )

            xml = re.sub(
                r"<[^>]+>",
                "",
                xml
            )

            xml = (
                xml.replace("&amp;", "&")
                   .replace("&lt;", "<")
                   .replace("&gt;", ">")
                   .replace("&quot;", '"')
                   .replace("&apos;", "'")
            )

            parts.append(xml)

    return "\n".join(parts)


def read_text(path):

    if path.suffix.lower() == ".docx":
        return read_docx(path)

    return path.read_text(
        encoding="utf-8",
        errors="replace"
    )


def words(text):

    return re.findall(
        r"\b[\w’'-]+\b",
        text,
        flags=re.UNICODE
    )


def canon_heading(x):

    x = x.strip().upper()

    x = re.sub(
        r"\s+",
        " ",
        x
    )

    return x


class Audit:

    def __init__(self):
        self.rows = []

    def add(
        self,
        category,
        check_id,
        status,
        message,
        observed="",
        expected=""
    ):

        self.rows.append({
            "category": category,
            "check_id": check_id,
            "status": status,
            "message": message,
            "observed": str(observed),
            "expected": str(expected),
        })

    def has_fail(self):

        return any(
            x["status"] == "FAIL"
            for x in self.rows
        )


# ============================================================
# Section parsing
# ============================================================

MAIN_HEADINGS = [
    "ABSTRACT",
    "BACKGROUND & SUMMARY",
    "METHODS",
    "DATA RECORDS",
    "DATA OVERVIEW",
    "TECHNICAL VALIDATION",
    "USAGE NOTES",
    "DATA AVAILABILITY",
    "CODE AVAILABILITY",
    "REFERENCES",
    "AUTHOR CONTRIBUTIONS",
    "ACKNOWLEDGEMENTS",
    "ACKNOWLEDGMENTS",
    "FUNDING",
    "COMPETING INTERESTS",
    "ETHICS",
    "ETHICS STATEMENT",
    "FIGURE LEGENDS",
]


def heading_positions(text):

    lines = text.splitlines()

    pos = []

    offset = 0

    known = {
        canon_heading(x): x
        for x in MAIN_HEADINGS
    }

    for line in lines:

        stripped = line.strip()
        c = canon_heading(stripped)

        if c in known:
            pos.append(
                (
                    known[c],
                    offset
                )
            )

        offset += len(line) + 1

    return pos


def section_map(text):

    pos = heading_positions(text)

    sections = {}

    for i, (name, start) in enumerate(pos):

        line_end = text.find(
            "\n",
            start
        )

        if line_end == -1:
            line_end = start

        body_start = line_end + 1

        end = (
            pos[i + 1][1]
            if i + 1 < len(pos)
            else len(text)
        )

        sections[name] = text[
            body_start:end
        ].strip()

    return sections


# ============================================================
# Title and abstract
# ============================================================

def audit_title_abstract(a, text, sections):

    first_nonempty = next(
        (
            x.strip()
            for x in text.splitlines()
            if x.strip()
        ),
        ""
    )

    nchar = len(first_nonempty)

    a.add(
        "title_abstract",
        "title_length",
        "PASS"
        if nchar <= 110
        else "FAIL",
        "Title length",
        nchar,
        "<=110 characters including spaces"
    )

    bad_title = re.findall(
        r"\b(?:novel|first|unique|unprecedented|"
        r"comprehensive|high-resolution)\b",
        first_nonempty,
        flags=re.I
    )

    a.add(
        "title_abstract",
        "title_claims",
        "WARN"
        if bad_title
        else "PASS",
        "Promotional/priority wording in title",
        bad_title,
        "none preferred"
    )

    abstract = sections.get(
        "ABSTRACT",
        ""
    )

    if not abstract:

        a.add(
            "title_abstract",
            "abstract_present",
            "FAIL",
            "Abstract section missing",
            "",
            "present"
        )

        return

    a.add(
        "title_abstract",
        "abstract_present",
        "PASS",
        "Abstract section present",
        "",
        ""
    )

    nwords = len(
        words(abstract)
    )

    a.add(
        "title_abstract",
        "abstract_word_count",
        "PASS"
        if nwords <= 170
        else "FAIL",
        "Abstract word count",
        nwords,
        "<=170 words"
    )

    if re.search(
        r"https?://|www\.",
        abstract,
        flags=re.I
    ):

        a.add(
            "title_abstract",
            "abstract_urls",
            "WARN",
            "URL detected in Abstract",
            "present",
            "none preferred"
        )

    else:

        a.add(
            "title_abstract",
            "abstract_urls",
            "PASS",
            "No URL in Abstract",
            "",
            ""
        )


# ============================================================
# Core sections
# ============================================================

def audit_sections(a, sections):

    required = [
        "ABSTRACT",
        "BACKGROUND & SUMMARY",
        "METHODS",
        "DATA RECORDS",
        "TECHNICAL VALIDATION",
        "DATA AVAILABILITY",
        "CODE AVAILABILITY",
        "REFERENCES",
    ]

    for x in required:

        a.add(
            "sections",
            "required_" + re.sub(
                r"\W+",
                "_",
                x.lower()
            ).strip("_"),
            "PASS"
            if x in sections
            else "FAIL",
            f"Required core section: {x}",
            "present"
            if x in sections
            else "missing",
            "present"
        )

    optional_but_useful = [
        "USAGE NOTES",
    ]

    for x in optional_but_useful:

        a.add(
            "sections",
            "recommended_" + re.sub(
                r"\W+",
                "_",
                x.lower()
            ).strip("_"),
            "PASS"
            if x in sections
            else "WARN",
            f"Recommended section: {x}",
            "present"
            if x in sections
            else "missing",
            "present recommended"
        )

    administrative = {
        "AUTHOR CONTRIBUTIONS":
            "Author contribution statement",
        "COMPETING INTERESTS":
            "Competing interests statement",
        "FUNDING":
            "Funding statement",
        "ACKNOWLEDGEMENTS":
            "Acknowledgements statement",
    }

    for x, label in administrative.items():

        present = x in sections

        if (
            x == "ACKNOWLEDGEMENTS"
            and "ACKNOWLEDGMENTS" in sections
        ):
            present = True

        a.add(
            "submission_statements",
            re.sub(
                r"\W+",
                "_",
                x.lower()
            ).strip("_"),
            "PASS"
            if present
            else "WARN",
            label,
            "present"
            if present
            else "missing",
            "add before submission if applicable"
        )


# ============================================================
# Methods result leakage
# ============================================================

def audit_methods(a, sections):

    s = sections.get(
        "METHODS",
        ""
    )

    if not s:
        return

    # Strong result-like frozen outcomes that are better placed
    # in Technical Validation / Data Records rather than Methods.
    result_numbers = {
        "15322": r"\b15[,]?322\b",
        "9539": r"\b9[,]?539\b",
        "476": r"\b476\b",
    }

    found = []

    for label, pat in result_numbers.items():

        if re.search(
            pat,
            s
        ):
            found.append(label)

    a.add(
        "methods",
        "result_statistics_in_methods",
        "WARN"
        if found
        else "PASS",
        (
            "Result-like final gene counts inside Methods"
        ),
        found,
        (
            "prefer procedures/thresholds in Methods; "
            "final counts in Technical Validation"
        )
    )

    result_phrases = [
        r"\bidentified\s+\d",
        r"\bshowed\s+(?:a|an|\d)",
        r"\bexplained\s+\d",
        r"\bthe largest\b",
        r"\bthe maximum\b",
        r"\bthe minimum\b",
    ]

    snippets = []

    for pat in result_phrases:

        for m in re.finditer(
            pat,
            s,
            flags=re.I
        ):

            snippets.append(
                re.sub(
                    r"\s+",
                    " ",
                    s[
                        max(0, m.start()-80):
                        min(len(s), m.end()+120)
                    ]
                )
            )

    a.add(
        "methods",
        "result_language",
        "WARN"
        if snippets
        else "PASS",
        "Result-like language in Methods",
        " | ".join(
            snippets[:5]
        ),
        "Methods should emphasize procedures"
    )


# ============================================================
# Data Records
# ============================================================

def audit_data_records(a, sections):

    s = sections.get(
        "DATA RECORDS",
        ""
    )

    if not s:
        return

    result_like = {
        "LRT 15,322":
            r"\b15[,]?322\b",
        "high-confidence 9,539":
            r"\b9[,]?539\b",
        "476 specificity genes":
            r"\b476\b",
        "PERMANOVA":
            r"\bPERMANOVA\b",
        "PCA":
            r"\bPC[12]\b|\bprincipal-component\b",
        "correlation result":
            r"\bSpearman\s+[ρr]?\s*=",
    }

    hit = []

    for name, pat in result_like.items():

        if re.search(
            pat,
            s,
            flags=re.I
        ):
            hit.append(name)

    a.add(
        "data_records",
        "result_content",
        "WARN"
        if hit
        else "PASS",
        "Result-like analytical statistics in Data Records",
        hit,
        (
            "Data Records should mainly describe files, "
            "formats and fields"
        )
    )

    # File/data terminology should be present.
    reusable_terms = [
        "metadata",
        "matrix",
        "Supplementary Table",
    ]

    missing = [
        x
        for x in reusable_terms
        if x.lower() not in s.lower()
    ]

    a.add(
        "data_records",
        "record_descriptions",
        "PASS"
        if not missing
        else "WARN",
        "Data-record descriptions present",
        f"missing={missing}",
        "metadata/matrices/tables described"
    )


# ============================================================
# Technical validation
# ============================================================

def audit_technical_validation(a, sections):

    s = sections.get(
        "TECHNICAL VALIDATION",
        ""
    )

    if not s:
        return

    concepts = {
        "sequencing_quality":
            r"\bQ30\b|\bsequencing quality\b",
        "alignment":
            r"\balignment\b|\bHISAT2\b",
        "gene_assignment":
            r"\bfeatureCounts\b|\bassignment rate\b",
        "gene_body_coverage":
            r"\bgene-body coverage\b",
        "expression_structure":
            r"\bPCA\b|\bprincipal-component\b",
        "permanova":
            r"\bPERMANOVA\b",
        "permdisp":
            r"\bPERMDISP\b",
        "technical_sensitivity":
            r"\bsensitivity\b|\bstress test\b",
        "sample_retention":
            r"\ball 77\b|\bretention\b",
        "historical_reference":
            r"\bhistorical-reference\b",
        "cross_cohort":
            r"\bcross-cohort\b",
    }

    missing = []

    for name, pat in concepts.items():

        if not re.search(
            pat,
            s,
            flags=re.I
        ):
            missing.append(name)

    a.add(
        "technical_validation",
        "coverage",
        "PASS"
        if not missing
        else "WARN",
        "Technical-validation coverage",
        f"missing={missing}",
        "major QC/robustness components represented"
    )


# ============================================================
# Usage Notes
# ============================================================

def audit_usage_notes(a, sections):

    s = sections.get(
        "USAGE NOTES",
        ""
    )

    if not s:
        return

    concepts = {
        "baseline_vs_0h":
            (
                r"\bBaseline\b.{0,120}\b0\s*h\b"
                r"|\b0\s*h\b.{0,120}\bBaseline\b"
            ),

        "historical_not_direct":
            (
                r"\bhistorical-reference\b.{0,180}"
                r"\b(?:not|differs)\b"
            ),

        "72h_caution":
            (
                r"\b72\s*h\b.{0,160}"
                r"\bhighest\b"
            ),

        "raw_counts_usage":
            (
                r"\braw gene-count\b"
                r"|\braw count\b"
            ),

        "TPM_usage":
            r"\bTPM\b",

        "VST_usage":
            r"\bVST\b|\bvariance-stabilized\b",
    }

    for name, pat in concepts.items():

        ok = bool(
            re.search(
                pat,
                s,
                flags=re.I | re.S
            )
        )

        a.add(
            "usage_notes",
            name,
            "PASS"
            if ok
            else "WARN",
            f"Usage guidance: {name}",
            "present"
            if ok
            else "missing",
            "present"
        )


# ============================================================
# References
# ============================================================

def audit_references(a, sections):

    s = sections.get(
        "REFERENCES",
        ""
    )

    if not s:
        return

    nums = [
        int(x)
        for x in re.findall(
            r"(?m)^\s*(\d+)\.\s+",
            s
        )
    ]

    if nums:

        expected = list(
            range(
                1,
                max(nums) + 1
            )
        )

        ok = nums == expected

    else:

        expected = []
        ok = False

    a.add(
        "references",
        "continuous_numbering",
        "PASS"
        if ok
        else "FAIL",
        "Reference numbering",
        nums,
        expected
    )

    dois = re.findall(
        r"10\.\d{4,9}/[-._;()/:A-Za-z0-9]+",
        s
    )

    a.add(
        "references",
        "doi_presence",
        "PASS"
        if len(dois) >= 8
        else "WARN",
        "DOIs provided for references where available",
        len(dois),
        "most journal references"
    )


# ============================================================
# Figure legends
# ============================================================

def audit_legends(a, sections):

    s = sections.get(
        "FIGURE LEGENDS",
        ""
    )

    if not s:

        a.add(
            "figure_legends",
            "section_present",
            "FAIL",
            "Figure legends section missing",
            "",
            "Figures 1-5"
        )

        return

    found = sorted({
        int(x)
        for x in re.findall(
            r"(?m)^Figure\s+([1-9]\d*)\s*\|",
            s,
            flags=re.I
        )
    })

    a.add(
        "figure_legends",
        "main_figures",
        "PASS"
        if found == [1, 2, 3, 4, 5]
        else "FAIL",
        "Main figure legends",
        found,
        [1, 2, 3, 4, 5]
    )


# ============================================================
# Placeholders
# ============================================================

def audit_placeholders(a, text):

    patterns = {
        "ethics":
            r"\[TO BE INSERTED\]",

        "geo":
            r"\bGSE[Xx]+\b",

        "github":
            r"\[repository URL[^\]]*\]",

        "zenodo":
            r"\[DOI to be inserted[^\]]*\]",

        "repository":
            r"\[repository/DOI[^\]]*\]",
    }

    for name, pat in patterns.items():

        n = len(
            re.findall(
                pat,
                text,
                flags=re.I
            )
        )

        a.add(
            "placeholders",
            name,
            "WARN"
            if n
            else "PASS",
            f"Unresolved {name} placeholder",
            n,
            0
        )


# ============================================================
# Output
# ============================================================

def write_outputs(a, outdir):

    outdir.mkdir(
        parents=True,
        exist_ok=True
    )

    detail = (
        outdir
        / "17B_scientific_data_structure_report.tsv"
    )

    with detail.open(
        "w",
        encoding="utf-8",
        newline=""
    ) as f:

        w = csv.DictWriter(
            f,
            delimiter="\t",
            fieldnames=[
                "category",
                "check_id",
                "status",
                "message",
                "observed",
                "expected",
            ],
            lineterminator="\n"
        )

        w.writeheader()
        w.writerows(a.rows)

    counts = defaultdict(
        lambda: {
            "PASS": 0,
            "WARN": 0,
            "FAIL": 0
        }
    )

    for r in a.rows:
        counts[
            r["category"]
        ][
            r["status"]
        ] += 1

    summary = (
        outdir
        / "17B_scientific_data_structure_summary.tsv"
    )

    with summary.open(
        "w",
        encoding="utf-8",
        newline=""
    ) as f:

        w = csv.writer(
            f,
            delimiter="\t",
            lineterminator="\n"
        )

        w.writerow([
            "category",
            "pass",
            "warn",
            "fail",
            "status"
        ])

        for cat in sorted(counts):

            c = counts[cat]

            status = (
                "FAIL"
                if c["FAIL"]
                else (
                    "WARN"
                    if c["WARN"]
                    else "PASS"
                )
            )

            w.writerow([
                cat,
                c["PASS"],
                c["WARN"],
                c["FAIL"],
                status
            ])

    P = sum(
        r["status"] == "PASS"
        for r in a.rows
    )

    W = sum(
        r["status"] == "WARN"
        for r in a.rows
    )

    F = sum(
        r["status"] == "FAIL"
        for r in a.rows
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

    log = (
        outdir
        / "17B_scientific_data_structure.log"
    )

    with log.open(
        "w",
        encoding="utf-8"
    ) as f:

        f.write("=" * 78 + "\n")
        f.write(
            "SCIENTIFIC DATA STRUCTURE & CONTENT AUDIT\n"
        )
        f.write("=" * 78 + "\n")

        f.write(
            f"OVERALL: {overall}\n"
        )

        f.write(
            f"PASS={P} WARN={W} FAIL={F}\n"
        )

        for r in a.rows:

            if r["status"] in {
                "WARN",
                "FAIL"
            }:

                f.write(
                    f"[{r['status']}] "
                    f"{r['category']} / "
                    f"{r['check_id']}: "
                    f"{r['message']} | "
                    f"observed={r['observed']} | "
                    f"expected={r['expected']}\n"
                )

    print("=" * 78)
    print(
        "SCIENTIFIC DATA STRUCTURE & CONTENT AUDIT"
    )
    print("=" * 78)

    print(
        f"OVERALL: {overall}   "
        f"PASS={P} WARN={W} FAIL={F}"
    )

    for cat in sorted(counts):

        c = counts[cat]

        status = (
            "FAIL"
            if c["FAIL"]
            else (
                "WARN"
                if c["WARN"]
                else "PASS"
            )
        )

        print(
            f"{cat:<25} "
            f"{status:<4} "
            f"PASS={c['PASS']:<3} "
            f"WARN={c['WARN']:<3} "
            f"FAIL={c['FAIL']:<3}"
        )

    print("-" * 78)
    print("Detailed:", detail)
    print("Summary :", summary)
    print("Log     :", log)
    print("=" * 78)

    return overall


# ============================================================
# Main
# ============================================================

def main():

    ap = argparse.ArgumentParser(
        description=(
            "Scientific Data Data Descriptor "
            "structure/content audit"
        )
    )

    ap.add_argument(
        "--manuscript",
        required=True
    )

    ap.add_argument(
        "--output-dir",
        default="17B_scientific_data_structure_audit"
    )

    args = ap.parse_args()

    p = Path(
        args.manuscript
    )

    if not p.is_file():

        print(
            "ERROR: manuscript not found:",
            p,
            file=sys.stderr
        )

        return 2

    try:

        text = normalize(
            read_text(p)
        )

    except Exception as e:

        print(
            "ERROR:",
            e,
            file=sys.stderr
        )

        return 2

    sections = section_map(
        text
    )

    a = Audit()

    a.add(
        "input",
        "manuscript_read",
        "PASS",
        "Manuscript loaded",
        f"{len(text):,} characters",
        ">1000 characters"
    )

    audit_title_abstract(
        a,
        text,
        sections
    )

    audit_sections(
        a,
        sections
    )

    audit_methods(
        a,
        sections
    )

    audit_data_records(
        a,
        sections
    )

    audit_technical_validation(
        a,
        sections
    )

    audit_usage_notes(
        a,
        sections
    )

    audit_references(
        a,
        sections
    )

    audit_legends(
        a,
        sections
    )

    audit_placeholders(
        a,
        text
    )

    overall = write_outputs(
        a,
        Path(
            args.output_dir
        )
    )

    return (
        1
        if overall == "FAIL"
        else 0
    )


if __name__ == "__main__":
    sys.exit(main())
