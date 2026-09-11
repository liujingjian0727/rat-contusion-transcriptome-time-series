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
# Frozen manuscript facts
# ============================================================

FROZEN = {
    "primary_samples": 77,
    "primary_baseline": 7,
    "primary_injury": 70,

    "historical_samples": 60,
    "historical_baseline": 6,
    "historical_injury": 54,

    "primary_filtered_genes": 18364,
    "historical_filtered_genes": 17328,
    "common_genes": 16918,

    "injury_only_filtered_genes": 18206,
    "lrt_dynamic_genes": 15322,
    "high_confidence_dynamic": 9539,
    "temporal_specificity_genes": 476,
}


EXPECTED_MAIN_FIGURES = [1, 2, 3, 4, 5]
EXPECTED_SUPP_FIGURES = list(range(1, 9))
EXPECTED_SUPP_TABLES = list(range(1, 9))


# ============================================================
# Phrases that should not survive into final manuscript
# ============================================================

FORBIDDEN_PATTERNS = {
    "validation_cohort": [
        r"\bvalidation cohort\b",
        r"\bexternal validation cohort\b",
        r"\bindependent validation cohort\b",
    ],

    "identical_cohort": [
        r"\bidentical cohort\b",
        r"\bidentical validation\b",
        r"\bexperimentally identical\b",
    ],

    "batch_correction_claim": [
        r"\bComBat[- ]corrected\b",
        r"\bbatch[- ]corrected cohorts\b",
        r"\bpooled cohorts\b",
        r"\bjointly normalized cohorts\b",
    ],

    "baseline_zero_equivalence": [
        r"\b0\s*h\s*(?:was|is|served as)\s+(?:the\s+)?(?:control|baseline)\b",
        r"\bBaseline\s*(?:=|was equivalent to)\s*0\s*h\b",
    ],

    "72h_overclaim": [
        r"\b72\s*h[- ]specific\b",
        r"\bspecific to 72\s*h\b",
        r"\bintrinsically specific to 72\s*h\b",
    ],

    "old_sample_number": [
        r"\b78 samples\b",
        r"\bn\s*=\s*78\b",
    ],

    "excluded_sample_present": [
        r"\bR2h_1\b",
    ],
}


# ============================================================
# Placeholders allowed before final deposition, but reported
# ============================================================

PLACEHOLDER_PATTERNS = {
    "ethics_placeholder": [
        r"\[TO BE INSERTED\]",
        r"ethics[^.\n]{0,80}\[[^\]]+\]",
    ],

    "geo_placeholder": [
        r"\bGSE[Xx]+\b",
        r"\[GSE[Xx]+\]",
    ],

    "github_placeholder": [
        r"\[repository URL[^\]]*\]",
        r"\bGitHub:\s*\[[^\]]+\]",
    ],

    "zenodo_placeholder": [
        r"\[DOI to be inserted[^\]]*\]",
        r"\bZenodo DOI:\s*\[[^\]]+\]",
    ],

    "generic_placeholder": [
        r"\bTBD\b",
        r"\bTODO\b",
        r"\bPLACEHOLDER\b",
        r"\[insert[^\]]*\]",
    ],
}


# ============================================================
# Required wording / concepts
# ============================================================

REQUIRED_CONCEPTS = {
    "historical_reference": [
        r"\bhistorical[- ]reference cohort\b",
        r"\bhistorical reference\b",
    ],

    "baseline_uninjured": [
        r"\bBaseline\b.{0,80}\buninjured\b",
        r"\buninjured\b.{0,80}\bBaseline\b",
    ],

    "zero_immediate": [
        r"\b0\s*h\b.{0,100}\bimmediately\b",
        r"\bimmediate(?:ly)? post[- ]injury\b",
    ],

    "no_batch_correction": [
        r"\bno\b.{0,80}\bbatch correction\b",
        r"\bno\b.{0,80}\bComBat\b",
        r"\bnot\b.{0,80}\bbatch[- ]correct",
    ],

    "terminal_72h_caution": [
        r"\b72\s*h\b.{0,150}\bhighest observed\b",
        r"\bterminal\b.{0,80}\b72\s*h\b",
        r"\bhighest\b.{0,80}\bwithin the sampled interval\b",
    ],
}


# ============================================================
# Helpers
# ============================================================

def normalize_text(text):
    text = text.replace("\u00a0", " ")
    text = text.replace("–", "-")
    text = text.replace("—", "-")
    text = text.replace("−", "-")
    return text


def read_docx(path):
    """
    Read visible text from DOCX using only Python stdlib.
    """
    parts = []

    with zipfile.ZipFile(path, "r") as z:

        targets = [
            "word/document.xml",
            "word/footnotes.xml",
            "word/endnotes.xml",
        ]

        for name in targets:

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

            # Basic XML entities
            xml = (
                xml.replace("&amp;", "&")
                   .replace("&lt;", "<")
                   .replace("&gt;", ">")
                   .replace("&quot;", '"')
                   .replace("&apos;", "'")
            )

            parts.append(xml)

    return "\n".join(parts)


def read_manuscript(path):

    suffix = path.suffix.lower()

    if suffix == ".docx":
        return read_docx(path)

    if suffix in {
        ".txt",
        ".md",
        ".tex",
        ".rtf"
    }:
        return path.read_text(
            encoding="utf-8",
            errors="replace"
        )

    raise ValueError(
        f"Unsupported manuscript type: {suffix}. "
        "Use .docx, .txt, .md, .tex or .rtf."
    )


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

    def failed(self):
        return any(
            r["status"] == "FAIL"
            for r in self.rows
        )


# ============================================================
# Figure/table reference extraction
# ============================================================

def main_figure_refs(text):

    refs = set()

    patterns = [
        r"\bFig(?:ure)?\.?\s*([1-5])\b",
        r"\bFigures?\s*([1-5])\b",
    ]

    for pat in patterns:
        for m in re.finditer(
            pat,
            text,
            flags=re.I
        ):
            refs.add(int(m.group(1)))

    return sorted(refs)


def supplementary_figure_refs(text):

    refs = set()

    patterns = [
        r"\bSupplementary\s+Fig(?:ure)?\.?\s*S?(\d+)\b",
        r"\bSupplementary\s+Figures?\s*S?(\d+)\b",
        r"\bFig(?:ure)?\.?\s*S(\d+)\b",
    ]

    for pat in patterns:
        for m in re.finditer(
            pat,
            text,
            flags=re.I
        ):
            refs.add(int(m.group(1)))

    return sorted(refs)


def supplementary_table_refs(text):

    refs = set()

    patterns = [
        r"\bSupplementary\s+Table\s*S?(\d+)\b",
        r"\bSupplementary\s+Tables\s*S?(\d+)\b",
        r"\bTable\s*S(\d+)\b",
    ]

    for pat in patterns:
        for m in re.finditer(
            pat,
            text,
            flags=re.I
        ):
            refs.add(int(m.group(1)))

    return sorted(refs)


# ============================================================
# Frozen-number checks
# ============================================================

def count_number(text, n):

    # Accept 18,364 and 18364.
    raw = str(n)
    comma = f"{n:,}"

    pat = (
        rf"(?<!\d)(?:"
        rf"{re.escape(raw)}|"
        rf"{re.escape(comma)}"
        rf")(?!\d)"
    )

    return len(
        re.findall(
            pat,
            text
        )
    )


def check_frozen_numbers(a, text):

    required = [
        ("primary_samples", 77),
        ("historical_samples", 60),
        ("primary_filtered_genes", 18364),
        ("historical_filtered_genes", 17328),
        ("common_genes", 16918),
        ("lrt_dynamic_genes", 15322),
        ("high_confidence_dynamic", 9539),
        ("temporal_specificity_genes", 476),
    ]

    for name, n in required:

        k = count_number(
            text,
            n
        )

        a.add(
            "frozen_numbers",
            name,
            "PASS" if k > 0 else "FAIL",
            f"Frozen number {n:,} appears in manuscript",
            k,
            ">=1 occurrence"
        )


# ============================================================
# Specific stage sample-size checks
# ============================================================

def check_primary_design(a, text):

    expected = {
        "Baseline": 7,
        "0": 9,
        "1": 8,
        "2": 7,
        "6": 6,
        "10": 6,
        "14": 7,
        "18": 6,
        "36": 7,
        "60": 7,
        "72": 7,
    }

    # We do not require every n to be repeated many times;
    # this checks whether the full design statement can be located.
    found = 0
    missing = []

    for stage, n in expected.items():

        if stage == "Baseline":

            pat = (
                rf"\bBaseline\b"
                rf"[^.\n]{{0,80}}"
                rf"(?:n\s*=\s*)?{n}\b"
            )

        else:

            pat = (
                rf"\b{stage}\s*h\b"
                rf"[^.\n]{{0,50}}"
                rf"(?:n\s*=\s*)?{n}\b"
            )

        if re.search(
            pat,
            text,
            flags=re.I
        ):
            found += 1

        else:
            missing.append(
                f"{stage}:{n}"
            )

    # Design may only be enumerated once.
    a.add(
        "study_design",
        "primary_stage_n",
        "PASS" if not missing else "WARN",
        "2026 stage-specific sample sizes",
        f"matched={found}/11; missing={missing}",
        "Baseline7 + 70 injury samples"
    )


def check_historical_design(a, text):

    expected = {
        "Baseline": 6,
        "4": 7,
        "8": 8,
        "12": 8,
        "16": 8,
        "20": 7,
        "24": 8,
        "48": 8,
    }

    found = 0
    missing = []

    for stage, n in expected.items():

        if stage == "Baseline":

            pat = (
                rf"\bBaseline\b"
                rf"[^.\n]{{0,80}}"
                rf"(?:n\s*=\s*)?{n}\b"
            )

        else:

            pat = (
                rf"\b{stage}\s*h\b"
                rf"[^.\n]{{0,50}}"
                rf"(?:n\s*=\s*)?{n}\b"
            )

        if re.search(
            pat,
            text,
            flags=re.I
        ):
            found += 1
        else:
            missing.append(
                f"{stage}:{n}"
            )

    a.add(
        "study_design",
        "historical_stage_n",
        "PASS" if not missing else "WARN",
        "2021 stage-specific sample sizes",
        f"matched={found}/8; missing={missing}",
        "60-sample historical-reference design"
    )


# ============================================================
# Forbidden phrases
# ============================================================

def check_forbidden(a, text):

    for group, patterns in FORBIDDEN_PATTERNS.items():

        matches = []

        for pat in patterns:

            for m in re.finditer(
                pat,
                text,
                flags=re.I
            ):

                left = max(
                    0,
                    m.start() - 70
                )

                right = min(
                    len(text),
                    m.end() + 70
                )

                snippet = re.sub(
                    r"\s+",
                    " ",
                    text[left:right]
                )

                matches.append(
                    snippet
                )

        a.add(
            "terminology",
            group,
            "FAIL" if matches else "PASS",
            "Forbidden/legacy wording check",
            " | ".join(matches[:5]),
            "no occurrences"
        )


# ============================================================
# Required concepts
# ============================================================

def check_required_concepts(a, text):

    for group, patterns in REQUIRED_CONCEPTS.items():

        found = any(
            re.search(
                pat,
                text,
                flags=re.I | re.S
            )
            for pat in patterns
        )

        a.add(
            "terminology",
            group,
            "PASS" if found else "WARN",
            "Required manuscript concept",
            found,
            "present"
        )


# ============================================================
# Placeholder checks
# ============================================================

def check_placeholders(a, text):

    for group, patterns in PLACEHOLDER_PATTERNS.items():

        matches = []

        for pat in patterns:

            matches.extend(
                m.group(0)
                for m in re.finditer(
                    pat,
                    text,
                    flags=re.I
                )
            )

        a.add(
            "placeholders",
            group,
            "WARN" if matches else "PASS",
            "Unresolved placeholder check",
            "; ".join(matches[:10]),
            "none before final submission"
        )


# ============================================================
# Figure/table citation consistency
# ============================================================

def check_references(a, text):

    mf = main_figure_refs(text)
    sf = supplementary_figure_refs(text)
    st = supplementary_table_refs(text)

    def audit_set(
        category,
        check_id,
        observed,
        expected
    ):

        missing = sorted(
            set(expected) - set(observed)
        )

        extra = sorted(
            set(observed) - set(expected)
        )

        a.add(
            "cross_references",
            check_id,
            "PASS"
            if not missing and not extra
            else "FAIL",
            category,
            f"observed={observed}; missing={missing}; extra={extra}",
            expected
        )

    audit_set(
        "Main figure references",
        "main_figures",
        mf,
        EXPECTED_MAIN_FIGURES
    )

    audit_set(
        "Supplementary figure references",
        "supplementary_figures",
        sf,
        EXPECTED_SUPP_FIGURES
    )

    audit_set(
        "Supplementary table references",
        "supplementary_tables",
        st,
        EXPECTED_SUPP_TABLES
    )


# ============================================================
# Legacy sample / cohort numbers
# ============================================================

def check_old_numbers(a, text):

    checks = {
        "53_samples": r"\b53\s+samples\b",
        "78_samples": r"\b78\s+samples\b",
        "R2h_1": r"\bR2h_1\b",
        "1154_masigpro": r"\b1154\b",
        "14619_old_filter": r"\b14[,]?619\b",
    }

    for name, pat in checks.items():

        m = list(
            re.finditer(
                pat,
                text,
                flags=re.I
            )
        )

        a.add(
            "legacy_values",
            name,
            "FAIL" if m else "PASS",
            "Legacy value/sample check",
            len(m),
            0
        )


# ============================================================
# Output
# ============================================================

def write_report(a, out_dir):

    out_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    detail = (
        out_dir
        / "17A_manuscript_consistency_report.tsv"
    )

    with detail.open(
        "w",
        encoding="utf-8",
        newline=""
    ) as f:

        writer = csv.DictWriter(
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
            lineterminator="\n",
        )

        writer.writeheader()
        writer.writerows(a.rows)

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
        out_dir
        / "17A_manuscript_consistency_summary.tsv"
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

            x = counts[cat]

            status = (
                "FAIL"
                if x["FAIL"]
                else (
                    "WARN"
                    if x["WARN"]
                    else "PASS"
                )
            )

            w.writerow([
                cat,
                x["PASS"],
                x["WARN"],
                x["FAIL"],
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
        out_dir
        / "17A_manuscript_consistency.log"
    )

    with log.open(
        "w",
        encoding="utf-8"
    ) as f:

        f.write("=" * 78 + "\n")
        f.write(
            "SCIENTIFIC DATA MANUSCRIPT CONSISTENCY AUDIT\n"
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
        "SCIENTIFIC DATA MANUSCRIPT CONSISTENCY AUDIT"
    )
    print("=" * 78)
    print(
        f"OVERALL: {overall}   "
        f"PASS={P} WARN={W} FAIL={F}"
    )

    for cat in sorted(counts):

        x = counts[cat]

        status = (
            "FAIL"
            if x["FAIL"]
            else (
                "WARN"
                if x["WARN"]
                else "PASS"
            )
        )

        print(
            f"{cat:<24} "
            f"{status:<4} "
            f"PASS={x['PASS']:<3} "
            f"WARN={x['WARN']:<3} "
            f"FAIL={x['FAIL']:<3}"
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
            "Scientific Data manuscript vs frozen "
            "figure/table/data consistency audit"
        )
    )

    ap.add_argument(
        "--manuscript",
        required=True,
        help=(
            "Manuscript file: .docx/.txt/.md/.tex/.rtf"
        )
    )

    ap.add_argument(
        "--output-dir",
        default="17A_manuscript_consistency_audit"
    )

    ap.add_argument(
        "--strict-placeholders",
        action="store_true",
        help=(
            "Treat unresolved placeholders as FAIL "
            "instead of WARN."
        )
    )

    args = ap.parse_args()

    manuscript = Path(
        args.manuscript
    )

    if not manuscript.is_file():
        print(
            "ERROR: manuscript not found:",
            manuscript,
            file=sys.stderr
        )
        return 2

    try:
        text = read_manuscript(
            manuscript
        )

    except Exception as e:
        print(
            "ERROR reading manuscript:",
            e,
            file=sys.stderr
        )
        return 2

    text = normalize_text(
        text
    )

    if len(text.strip()) < 1000:
        print(
            "ERROR: extracted manuscript text appears too short.",
            file=sys.stderr
        )
        return 2

    a = Audit()

    a.add(
        "input",
        "manuscript_read",
        "PASS",
        "Manuscript text extracted",
        f"{len(text):,} characters",
        ">1,000 characters"
    )

    check_references(
        a,
        text
    )

    check_frozen_numbers(
        a,
        text
    )

    check_primary_design(
        a,
        text
    )

    check_historical_design(
        a,
        text
    )

    check_forbidden(
        a,
        text
    )

    check_required_concepts(
        a,
        text
    )

    check_placeholders(
        a,
        text
    )

    check_old_numbers(
        a,
        text
    )

    # Upgrade placeholder WARN -> FAIL when manuscript is truly final.
    if args.strict_placeholders:

        for r in a.rows:

            if (
                r["category"] == "placeholders"
                and r["status"] == "WARN"
            ):
                r["status"] = "FAIL"

    overall = write_report(
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
