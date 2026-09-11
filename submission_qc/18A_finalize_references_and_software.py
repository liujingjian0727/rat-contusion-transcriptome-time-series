#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import re
import sys


FILE = Path("Scientific_Data_manuscript_working.txt")


def die(msg):
    raise SystemExit("ERROR: " + msg)


if not FILE.is_file():
    die(f"not found: {FILE}")

text = FILE.read_text(encoding="utf-8")


# ============================================================
# Split manuscript into body / references / figure legends
# ============================================================

marker_refs = "\nREFERENCES\n"
marker_figs = "\nFIGURE LEGENDS\n"

if marker_refs not in text:
    die("REFERENCES heading not found")

if marker_figs not in text:
    die("FIGURE LEGENDS heading not found")

body, rest = text.split(marker_refs, 1)
refs_old, legends = rest.split(marker_figs, 1)


# ============================================================
# 1. Renumber existing in-text citations
#
# Existing:
# 8  HISAT2
# 9  featureCounts
# 10 DESeq2
# 11 RSeQC
# 12 Anderson
#
# New:
# 8  GRCr8
# 9  HISAT2
# 10 featureCounts
# 11 Rsubread
# 12 TMM
# 13 DESeq2
# 14 RSeQC
# 15 Anderson
# 16 vegan
# ============================================================

old_to_new = {
    8: 9,
    9: 10,
    10: 13,
    11: 14,
    12: 15,
}


def renumber_citation(match):

    inner = match.group(1)

    # Only process simple numeric citation blocks such as:
    # [8], [4,5], [4, 6]
    if not re.fullmatch(
        r"\s*\d+(?:\s*,\s*\d+)*\s*",
        inner
    ):
        return match.group(0)

    nums = [
        int(x)
        for x in re.findall(r"\d+", inner)
    ]

    nums = [
        old_to_new.get(x, x)
        for x in nums
    ]

    return "[" + ",".join(map(str, nums)) + "]"


body = re.sub(
    r"\[([0-9,\s]+)\]",
    renumber_citation,
    body
)


# ============================================================
# 2. GRCr8 + HISAT2 + featureCounts + Rsubread citations
# ============================================================

old = (
    "Quality-controlled reads were aligned to the rat GRCr8 "
    "reference genome using HISAT2 v2.2.1 [9]. Gene-level read "
    "assignment was performed using featureCounts implemented "
    "in Rsubread v2.24.0 [10]."
)

new = (
    "Quality-controlled reads were aligned to the rat GRCr8 "
    "reference genome [8] using HISAT2 v2.2.1 [9]. Gene-level "
    "read assignment was performed using featureCounts [10] "
    "implemented in Rsubread v2.24.0 [11]."
)

if old not in body:
    die("alignment/Rsubread target sentence not found")

body = body.replace(old, new, 1)


# ============================================================
# 3. Add TMM citation
# ============================================================

old = (
    "A corresponding processed TMM-scaled TPM matrix contained "
    "the same genes and samples and was retained as a "
    "complementary expression representation."
)

new = (
    "A corresponding processed TMM-scaled TPM matrix, using the "
    "trimmed mean of M-values normalization framework [12], "
    "contained the same genes and samples and was retained as a "
    "complementary expression representation."
)

if old not in body:
    die("TMM target sentence not found")

body = body.replace(old, new, 1)


# ============================================================
# 4. PERMANOVA + vegan citation
# ============================================================

old = (
    "Multivariate sampling-stage effects were evaluated using "
    "PERMANOVA with 9,999 permutations. Homogeneity of "
    "multivariate dispersion was assessed using betadisper with "
    "spatial-median centroids and bias adjustment, followed by "
    "permutation testing."
)

new = (
    "Multivariate sampling-stage effects were evaluated using "
    "PERMANOVA [15] with 9,999 permutations, and homogeneity of "
    "multivariate dispersion was assessed using betadisper with "
    "spatial-median centroids and bias adjustment, followed by "
    "permutation testing. These procedures were implemented "
    "using vegan v2.7-5 [16]."
)

if old not in body:
    die("PERMANOVA/vegan target paragraph not found")

body = body.replace(old, new, 1)


# ============================================================
# 5. Make verified R-package versions explicit
# ============================================================

old = (
    "The final computational environment included R v4.5.3 and "
    "Bioconductor v3.22. Major software used in the workflow "
    "included fastp v0.23.2, HISAT2 v2.2.1, Rsubread v2.24.0, "
    "DESeq2, RSeQC, vegan and ggplot2."
)

new = (
    "The final computational environment included R v4.5.3 and "
    "Bioconductor v3.22. Major software used in the workflow "
    "included fastp v0.23.2, HISAT2 v2.2.1, Rsubread v2.24.0, "
    "DESeq2 v1.50.2, RSeQC, vegan v2.7-5 and ggplot2 v4.0.3."
)

if old not in body:
    die("Code Availability software-version sentence not found")

body = body.replace(old, new, 1)


# ============================================================
# 6. Preserve currently verified References 1-7 exactly
# ============================================================

existing = {}

for line in refs_old.splitlines():

    line = line.strip()

    m = re.match(
        r"^(\d+)\.\s+(.+)$",
        line
    )

    if m:
        existing[int(m.group(1))] = line


for i in range(1, 8):

    if i not in existing:
        die(f"existing reference {i} not found")


refs_1_to_7 = [
    existing[i]
    for i in range(1, 8)
]


# ============================================================
# 7. Final References 8-16
# ============================================================

refs_8_to_16 = [

    (
        "8. Li, K. et al. Construction and evaluation of a new "
        "rat reference genome assembly, GRCr8, from long reads "
        "and long-range scaffolding. Genome Research 34, "
        "2081–2093 (2024). "
        "https://doi.org/10.1101/gr.279292.124."
    ),

    (
        "9. Kim, D., Paggi, J. M., Park, C., Bennett, C. & "
        "Salzberg, S. L. Graph-based genome alignment and "
        "genotyping with HISAT2 and HISAT-genotype. "
        "Nature Biotechnology 37, 907–915 (2019). "
        "https://doi.org/10.1038/s41587-019-0201-4."
    ),

    (
        "10. Liao, Y., Smyth, G. K. & Shi, W. featureCounts: "
        "an efficient general-purpose program for assigning "
        "sequence reads to genomic features. Bioinformatics 30, "
        "923–930 (2014). "
        "https://doi.org/10.1093/bioinformatics/btt656."
    ),

    (
        "11. Liao, Y., Smyth, G. K. & Shi, W. The R package "
        "Rsubread is easier, faster, cheaper and better for "
        "alignment and quantification of RNA sequencing reads. "
        "Nucleic Acids Research 47, e47 (2019). "
        "https://doi.org/10.1093/nar/gkz114."
    ),

    (
        "12. Robinson, M. D. & Oshlack, A. A scaling "
        "normalization method for differential expression "
        "analysis of RNA-seq data. Genome Biology 11, R25 "
        "(2010). "
        "https://doi.org/10.1186/gb-2010-11-3-r25."
    ),

    (
        "13. Love, M. I., Huber, W. & Anders, S. Moderated "
        "estimation of fold change and dispersion for RNA-seq "
        "data with DESeq2. Genome Biology 15, 550 (2014). "
        "https://doi.org/10.1186/s13059-014-0550-8."
    ),

    (
        "14. Wang, L., Wang, S. & Li, W. RSeQC: quality "
        "control of RNA-seq experiments. Bioinformatics 28, "
        "2184–2185 (2012). "
        "https://doi.org/10.1093/bioinformatics/bts356."
    ),

    (
        "15. Anderson, M. J. A new method for non-parametric "
        "multivariate analysis of variance. Austral Ecology 26, "
        "32–46 (2001). "
        "https://doi.org/10.1111/"
        "j.1442-9993.2001.01070.pp.x."
    ),

    (
        "16. Oksanen, J. et al. vegan: Community Ecology "
        "Package. R package version 2.7-5 (2026). "
        "https://doi.org/10.32614/CRAN.package.vegan."
    ),
]


refs_final = refs_1_to_7 + refs_8_to_16


# ============================================================
# 8. Internal validation of references
# ============================================================

numbers = []

for x in refs_final:

    m = re.match(
        r"^(\d+)\.",
        x
    )

    if not m:
        die("malformed reference: " + x)

    numbers.append(
        int(m.group(1))
    )


if numbers != list(range(1, 17)):
    die(
        f"reference numbering invalid: {numbers}"
    )


# ============================================================
# 9. Verify citations used in body
# ============================================================

citation_numbers = set()

for m in re.finditer(
    r"\[([0-9,\s]+)\]",
    body
):

    for x in re.findall(
        r"\d+",
        m.group(1)
    ):
        citation_numbers.add(
            int(x)
        )


missing_citations = [
    x
    for x in range(1, 17)
    if x not in citation_numbers
]


if missing_citations:
    print(
        "WARNING: references not detected in body:",
        missing_citations
    )


invalid_citations = sorted(
    x
    for x in citation_numbers
    if x > 16
)


if invalid_citations:
    die(
        f"undefined citations detected: {invalid_citations}"
    )


# ============================================================
# 10. Reassemble
# ============================================================

new_text = (
    body.rstrip()
    + "\n\nREFERENCES\n\n"
    + "\n\n".join(refs_final)
    + "\n\nFIGURE LEGENDS\n\n"
    + legends.lstrip()
)

FILE.write_text(
    new_text,
    encoding="utf-8"
)


print("=" * 78)
print("STEP 18A FINAL REFERENCE + SOFTWARE PATCH")
print("=" * 78)
print("PASS: GRCr8 reference added")
print("PASS: HISAT2 renumbered to Ref. 9")
print("PASS: featureCounts renumbered to Ref. 10")
print("PASS: Rsubread Ref. 11 added")
print("PASS: TMM Ref. 12 added")
print("PASS: DESeq2 renumbered to Ref. 13")
print("PASS: RSeQC renumbered to Ref. 14")
print("PASS: Anderson PERMANOVA renumbered to Ref. 15")
print("PASS: vegan v2.7-5 added as Ref. 16")
print("PASS: software versions updated")
print("References: 1–16 continuous")
print(
    "Body citation numbers:",
    sorted(citation_numbers)
)

if missing_citations:
    print(
        "REVIEW: uncited references:",
        missing_citations
    )
else:
    print("PASS: all References 1–16 cited in manuscript")

print("Output:", FILE)
print("=" * 78)
