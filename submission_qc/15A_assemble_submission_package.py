#!/usr/bin/env python3

from pathlib import Path
import shutil
import hashlib
import csv
import re
import subprocess
import sys

# ============================================================
# 15A
# Assemble Scientific Data submission figure / table package
# Main Figures 1-5
# Supplementary Figures S1-S8
# Supplementary Tables S1-S8
# Figure Source Data
# Captions / legends / manifest / SHA256
# ============================================================

ROOT = Path.cwd()

PACKAGE = ROOT / "submission_package_scientific_data"

MAIN_DIR = PACKAGE / "01_Main_Figures"
SUPP_FIG_DIR = PACKAGE / "02_Supplementary_Figures"
SUPP_TABLE_DIR = PACKAGE / "03_Supplementary_Tables"
SOURCE_DIR = PACKAGE / "04_Source_Data"
CAPTION_DIR = PACKAGE / "05_Captions_and_Indexes"
MANIFEST_DIR = PACKAGE / "06_Manifests"

FINAL_FIG_DIR = ROOT / "final_submission_figures"
SUPP_ROOT = ROOT / "supplementary_figures"

CLEAN_EXISTING_PACKAGE = True


# ============================================================
# 0. Helpers
# ============================================================

def banner(text):
    print()
    print("=" * 60)
    print(text)
    print("=" * 60)


def sha256_file(path, chunk_size=1024 * 1024):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            chunk = f.read(chunk_size)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def copy_file(src, dst):
    src = Path(src)
    dst = Path(dst)

    if not src.exists():
        raise FileNotFoundError(str(src))

    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)

    return dst


def copy_optional(src, dst):
    src = Path(src)
    dst = Path(dst)

    if not src.exists():
        print(f"OPTIONAL MISSING: {src}")
        return None

    return copy_file(src, dst)


def safe_name(x):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", x)


def image_info(path):
    """
    Return ImageMagick identify summary when available.
    PDFs are skipped.
    """
    path = Path(path)

    if path.suffix.lower() not in {
        ".tif",
        ".tiff",
        ".png",
        ".jpg",
        ".jpeg"
    }:
        return ""

    try:
        result = subprocess.run(
            [
                "identify",
                "-format",
                "%wx%h; %x x %y; %[units]",
                str(path)
            ],
            capture_output=True,
            text=True,
            check=True
        )
        return result.stdout.strip()
    except Exception:
        return ""


def locate_unique_by_globs(base_dir, patterns, description):
    """
    Try glob patterns in order.
    A pattern is accepted only when it resolves to exactly one file.
    If a pattern gives >1 file, continue to more specific patterns only
    if earlier patterns were broad. At the end ambiguity is an error.
    """

    base_dir = Path(base_dir)

    all_hits = []

    for pattern in patterns:
        hits = sorted([
            x for x in base_dir.glob(pattern)
            if x.is_file()
        ])

        if len(hits) == 1:
            return hits[0]

        if len(hits) > 1:
            all_hits.extend(hits)

    unique_hits = sorted(set(all_hits))

    if len(unique_hits) == 1:
        return unique_hits[0]

    if len(unique_hits) > 1:
        print()
        print(f"AMBIGUOUS {description}:")
        for x in unique_hits:
            print(f"  {x}")
        raise RuntimeError(
            f"Multiple candidate files found for {description}."
        )

    raise FileNotFoundError(
        f"No file found for {description}."
    )


def recursive_find_unique(basenames=None, regex=None, required=True):
    """
    Locate one file recursively beneath ROOT.
    """

    files = [
        x for x in ROOT.rglob("*")
        if x.is_file()
        and PACKAGE not in x.parents
    ]

    if basenames:
        for name in basenames:
            hits = [
                x for x in files
                if x.name == name
            ]

            if len(hits) == 1:
                return hits[0]

            if len(hits) > 1:
                # Prefer current analysis tree / shallower path
                hits = sorted(
                    hits,
                    key=lambda p: (
                        len(p.parts),
                        str(p)
                    )
                )

                print(
                    f"NOTE: multiple '{name}' files found; "
                    f"using: {hits[0]}"
                )
                return hits[0]

    if regex:
        rx = re.compile(regex, flags=re.I)

        hits = [
            x for x in files
            if rx.search(x.name)
        ]

        if len(hits) == 1:
            return hits[0]

        if len(hits) > 1:
            hits = sorted(
                hits,
                key=lambda p: (
                    len(p.parts),
                    str(p)
                )
            )

            print(
                f"NOTE: multiple regex matches; using: {hits[0]}"
            )
            return hits[0]

    if required:
        raise FileNotFoundError(
            f"Could not locate required file. "
            f"basenames={basenames}, regex={regex}"
        )

    return None


def read_tsv(path):
    with open(path, "r", newline="") as f:
        return list(
            csv.DictReader(
                f,
                delimiter="\t"
            )
        )


def write_tsv(path, rows, fieldnames):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


# ============================================================
# 1. Prepare package
# ============================================================

banner(
    "15A ASSEMBLE SCIENTIFIC DATA SUBMISSION PACKAGE"
)

if CLEAN_EXISTING_PACKAGE and PACKAGE.exists():
    shutil.rmtree(PACKAGE)

for d in [
    MAIN_DIR,
    SUPP_FIG_DIR,
    SUPP_TABLE_DIR,
    SOURCE_DIR,
    CAPTION_DIR,
    MANIFEST_DIR
]:
    d.mkdir(
        parents=True,
        exist_ok=True
    )

print(f"Working directory : {ROOT}")
print(f"Package directory : {PACKAGE}")


# ============================================================
# 2. Main Figure captions
# ============================================================

MAIN_CAPTIONS = {

1: (
"Figure 1 | Experimental design and data-resource overview. "
"a, Overview of the 2026 primary rat skeletal-muscle contusion cohort, "
"comprising 7 uninjured Baseline samples and 70 post-injury samples "
"collected after a 70-cm gravity-hammer contusion, for a total of 77 "
"samples. b, Sampling-stage design spanning Baseline and ten post-injury "
"stages from 0 to 72 h. Sampling-stage positions are shown at equal "
"spacing for readability and do not represent proportional elapsed time. "
"c, Analysis framework connecting raw-count filtering and variance-"
"stabilized expression analyses with TMM-normalized TPM-based descriptive "
"temporal annotations and independent Baseline-relative cross-cohort "
"comparison. d, Cohort-level comparison between the 2026 primary resource "
"(70-cm injury, 77 samples, 0–72 h) and the previously published 2021 "
"historical-reference cohort (50-cm injury, 60 samples, 4–48 h). The two "
"cohorts were analysed independently and were not pooled or batch-corrected."
),

2: (
"Figure 2 | Internal technical validation of the 2026 primary "
"transcriptomic resource. a, Relationship between gene-assigned library "
"size and the number of detected genes across the 77 RNA-seq samples, "
"calculated from the unfiltered raw count matrix. Colours indicate sampling "
"stages. b, Principal component analysis of variance-stabilized expression "
"profiles for 18,364 genes retained after predefined count filtering. PC1 "
"and PC2 explained 28.28% and 12.85% of the total variance, respectively. "
"c, Pairwise Spearman correlation matrix of variance-stabilized expression "
"profiles across all 77 samples, ordered by sampling stage. White "
"boundaries delineate sampling-stage groups. d, Quantitative sample-level "
"expression diagnostics. Sixty-five samples were not flagged, seven were "
"identified for single-metric review, and five for multi-metric review. "
"Review status was used for additional sensitivity assessment rather than "
"automatic sample exclusion; all 77 samples were retained in the primary "
"analyses."
),

3: (
"Figure 3 | Temporal organization of the 2026 primary transcriptomic "
"resource. a, Principal component representation of variance-stabilized "
"expression profiles with chronological sampling-stage centroids connected "
"to summarize the global temporal trajectory. b, Global multivariate "
"separation among sampling stages. Sampling stage explained 41.5% of the "
"full variance-stabilized Euclidean expression space (PERMANOVA, R²=0.415, "
"P<0.001), whereas within-stage multivariate dispersion did not differ "
"significantly (PERMDISP, P=0.722). c, Root-mean-square distances between "
"consecutive stage centroids in the full 18,364-gene variance-stabilized "
"expression space. These values are descriptive measures of transcriptomic "
"transition magnitude and are not interpreted as biological velocities. "
"d, Euclidean distance of individual samples from their corresponding "
"sampling-stage centroid, providing a complementary view of within-stage "
"expression dispersion."
),

4: (
"Figure 4 | Differential and continuous temporal expression dynamics in "
"the 2026 primary cohort. a, Numbers of up- and downregulated genes at each "
"post-injury stage relative to the uninjured Baseline group, using DESeq2 "
"with adjusted P<0.05 and |log2 fold change|≥1. b, Hierarchical definition "
"of temporal gene sets: 18,364 genes after predefined count filtering; "
"18,206 genes retained by the injury-sample expression filter; 15,322 genes "
"significant in the continuous-time spline likelihood-ratio test at FDR<0.05; "
"9,539 genes additionally satisfying a mean variance-stabilized temporal "
"range of at least 1; and 476 high-confidence dynamic genes additionally "
"meeting injury-stage Tau≥0.85 and peak SPM≥0.50 criteria. c, Temporal "
"heatmap of the 476 genes using sampling-stage median log2(TMM-normalized "
"TPM+1) values, row-standardized for visualization and ordered by peak stage "
"and specificity. d, Distribution of the 476 genes by the stage at which "
"their highest observed injury-stage expression occurred. The 72-h stage "
"is the terminal sampled time point; genes peaking at 72 h are therefore "
"described as having their highest observed expression at 72 h rather than "
"as intrinsically 72-h-specific."
),

5: (
"Figure 5 | Cross-cohort concordance of Baseline-relative transcriptomic "
"responses. a, Spearman correlations between genome-wide Baseline-relative "
"log2 fold-change vectors for all 10 × 7 combinations of 2026 primary-"
"cohort and 2021 historical-reference stages, calculated across 16,918 "
"genes retained after independent count filtering in both cohorts. "
"b, Chronologically nearest cross-cohort stage comparisons; ties in absolute "
"sampling-time difference are retained. c, Direction concordance among "
"shared differentially expressed genes, with DEGs defined independently "
"within each cohort using adjusted P<0.05 and |log2 fold change|≥1. "
"d, Sensitivity of the cross-cohort concordance structure to removal of "
"reviewed samples. Each cohort was analysed relative to its own uninjured "
"Baseline group; no pooling or cross-cohort batch correction was performed. "
"Unconstrained high-concordance stage pairs are interpreted descriptively "
"and not as evidence of chronological or biological equivalence."
)
}


# ============================================================
# 3. Supplementary Figure captions
# ============================================================

SUPP_CAPTIONS = {

1: (
"Supplementary Figure S1 | Gene-body coverage profiles in the 2026 primary "
"cohort. a, RSeQC-derived gene-body coverage profiles for all 77 final "
"samples after exclusion of the non-final R2h_1 sample. Thin lines represent "
"individual samples, the black line represents the cohort median and the "
"shaded region represents the interquartile range. Profiles were normalized "
"within sample to the RSeQC-style 0–1 relative-coverage scale. b, Median "
"gene-body coverage profile for each sampling stage. Gene-body position is "
"shown from the 5′ to the 3′ end. These profiles were used as descriptive "
"technical-quality diagnostics and motivated explicit coverage-sensitivity "
"analyses rather than post-hoc sample exclusion."
),

2: (
"Supplementary Figure S2 | Sensitivity of global expression structure to "
"gene-body coverage variation. a, Principal component analysis of the "
"original variance-stabilized expression matrix. b, PCA after removal of "
"coverage-associated expression components conditional on sampling stage, "
"thereby preserving stage-associated effects. c, PCA after a more aggressive "
"coverage-only adjustment that removes all linear signal associated with the "
"coverage metrics without explicitly protecting sampling stage. d, Sampling-"
"stage PERMANOVA effect sizes in the original, stage-preserving and coverage-"
"only expression spaces. The corresponding R² values were 0.415, 0.417 and "
"0.334 (all P<0.001), with nonsignificant PERMDISP tests (P=0.722, 0.821 "
"and 0.894, respectively). The coverage-only analysis is interpreted as a "
"conservative sensitivity bound rather than a preferred expression matrix."
),

3: (
"Supplementary Figure S3 | Raw-read sequencing depth, alignment and gene-"
"level assignment quality control for the 2026 primary cohort. a, Number "
"of paired-end raw read pairs supplied by the sequencing provider across "
"the 77 final samples. b, HISAT2 overall alignment rates. c, Gene-level "
"assignment rates reported by featureCounts. d, Relationship between HISAT2 "
"overall alignment rate and featureCounts assignment rate. Individual points "
"represent samples and boxplots summarize sampling-stage distributions. "
"Sequencing depth and alignment efficiencies were highly consistent across "
"the cohort, whereas gene-level assignment rates showed modest variation "
"among sampling stages. No association was observed between alignment and "
"assignment efficiencies (Spearman ρ=-0.060, P=0.602)."
),

4: (
"Supplementary Figure S4 | Sensitivity of global transcriptomic structure "
"to fastp-derived sequencing-quality variation. a, Principal component "
"analysis of the original variance-stabilized expression matrix for the "
"77-sample 2026 primary cohort. b, PCA after removal of expression components "
"associated with Q30 base-call quality, GC content and estimated read "
"duplication conditional on sampling stage. This stage-preserving adjustment "
"retains stage-associated expression effects while removing residual "
"QC-associated variation. c, PCA following a more conservative QC-only "
"adjustment in which all linear expression components associated with Q30, "
"GC content and duplication were removed without explicitly preserving "
"sampling stage. Because these QC metrics were themselves associated with "
"sampling stage, the QC-only adjustment represents a sensitivity bound and "
"may also remove genuine stage-related biological variation. d, PERMANOVA "
"effect sizes for sampling stage in the original and adjusted full-expression "
"spaces. Sampling-stage structure remained significant in all three spaces "
"(all P<0.001), with R² values of 0.415, 0.539 and 0.273 for the original, "
"stage-preserving and QC-only spaces, respectively. PERMDISP remained "
"nonsignificant (P=0.722, 0.819 and 0.521). PCA was recalculated independently "
"for each matrix; principal-component orientation and sign are therefore not "
"directly comparable among panels."
),

5: (
"Supplementary Figure S5 | Internal transcriptomic structure of the "
"independently reprocessed 2021 historical-reference cohort. a, Principal "
"component analysis of variance-stabilized expression profiles across 60 "
"samples and 17,328 genes retained after predefined count filtering. Points "
"represent individual samples and arrows connect sampling-stage centroids in "
"chronological order. b, Pairwise Spearman correlation matrix of variance-"
"stabilized expression profiles across all 60 samples, ordered by sampling "
"stage. White boundaries delineate sampling-stage groups. c, Adjacent-stage "
"transcriptomic transition magnitudes quantified as root-mean-square "
"distances between consecutive stage centroids in the full variance-"
"stabilized expression space. These distances are descriptive rather than "
"biological velocities. d, Euclidean distance of individual samples from "
"their corresponding stage centroid. Global transcriptomic profiles differed "
"significantly among stages (PERMANOVA, R²=0.520, P<0.001), whereas "
"multivariate dispersion did not differ significantly (PERMDISP, P=0.237). "
"The 2021 dataset was analysed independently as a previously published "
"historical-reference cohort and was not pooled or batch-corrected with the "
"2026 primary cohort."
),

6: (
"Supplementary Figure S6 | Baseline-relative differential-expression "
"responses in the independently reprocessed 2021 historical-reference "
"cohort. a, Numbers of up- and downregulated genes at each post-injury stage "
"relative to the uninjured Baseline group. DEGs were defined using DESeq2 "
"with adjusted P<0.05 and |log2 fold change|≥1. b, Genome-wide magnitude of "
"Baseline-relative responses across 17,328 genes. Circles indicate median "
"absolute log2 fold change, error bars the interquartile range and triangles "
"the 90th percentile. c, Pairwise Spearman correlations between genome-wide "
"Baseline-relative log2 fold-change vectors for the seven post-injury stages. "
"d, Spearman correlations between genome-wide response vectors for consecutive "
"sampled stages. These analyses describe internal temporal organization of "
"the historical-reference cohort and are not interpreted as evidence of "
"temporal equivalence with the 2026 primary cohort."
),

7: (
"Supplementary Figure S7 | Complete cross-cohort diagnostics of Baseline-"
"relative transcriptomic responses in the 2026 primary cohort and the "
"independently reprocessed 2021 historical-reference cohort. a, Pearson "
"correlations between genome-wide Baseline-relative log2 fold-change vectors "
"for all 70 combinations of the ten 2026 and seven 2021 post-injury stages, "
"calculated across 16,918 common genes. b, Numbers of shared differentially "
"expressed genes for each cross-cohort stage pair. DEGs were defined "
"independently within each cohort using adjusted P<0.05 and |log2 fold change|"
"≥1. c, Direction concordance among shared DEGs, defined as the proportion "
"showing the same log2 fold-change sign in both cohorts. d, Comparison of "
"Spearman and Pearson genome-wide response concordance across all 70 stage "
"pairs; colours denote absolute sampling-time difference. The two cohorts "
"were analysed independently relative to their own Baseline groups and were "
"not pooled or batch-corrected. High concordance between chronologically "
"distant stages is descriptive and does not establish temporal equivalence."
),

8: (
"Supplementary Figure S8 | Robustness of temporal and cross-cohort "
"transcriptomic structure to exclusion of reviewed samples. a, Sampling-stage "
"PERMANOVA effect sizes in the complete 77-sample primary cohort and after "
"exclusion of five multi-metric expression-review samples, exclusion of the "
"sequencing-review sample R2h_4, or exclusion of the union of all six reviewed "
"samples. PERMDISP P values are shown for each scenario. b, Stability of the "
"complete 10 × 7 cross-cohort Spearman response-concordance matrix after "
"multi-sample exclusion, quantified by correlation with the primary matrix "
"and mean absolute matrix difference. c, Sampling-stage PERMANOVA effect "
"sizes after leave-one-out removal of each reviewed sample individually. "
"d, Correlation between each leave-one-out cross-cohort response-concordance "
"matrix and the primary 70-cell matrix; labels additionally report the mean "
"absolute matrix difference. All leave-one-out analyses retained significant "
"sampling-stage structure and no individual reviewed sample dominated the "
"cross-cohort response pattern. Multi-sample exclusion changed some "
"quantitative estimates and exact DEG counts and is therefore interpreted "
"as a sensitivity analysis rather than evidence that sample removal had no "
"effect."
)
}


# ============================================================
# 4. Locate main figures
# ============================================================

banner("LOCATING MAIN FIGURES")

if not FINAL_FIG_DIR.exists():
    raise FileNotFoundError(
        f"Missing directory: {FINAL_FIG_DIR}"
    )

MAIN_EXACT = {
    1: {
        "pdf": "Figure_1_experimental_design_and_resource_overview.pdf",
        "tiff": "Figure_1_experimental_design_and_resource_overview_600dpi.tiff"
    },
    2: {
        "pdf": "Figure_2_internal_technical_validation.pdf",
        "tiff": "Figure_2_internal_technical_validation_600dpi.tiff"
    }
}

main_sources = {}

for fig_no in range(1, 6):

    if fig_no in MAIN_EXACT:
        pdf = FINAL_FIG_DIR / MAIN_EXACT[fig_no]["pdf"]
        tif = FINAL_FIG_DIR / MAIN_EXACT[fig_no]["tiff"]

        if not pdf.exists():
            pdf = locate_unique_by_globs(
                FINAL_FIG_DIR,
                [
                    f"Figure_{fig_no}_*.pdf",
                    f"Figure{fig_no}_*.pdf"
                ],
                f"Figure {fig_no} PDF"
            )

        if not tif.exists():
            tif = locate_unique_by_globs(
                FINAL_FIG_DIR,
                [
                    f"Figure_{fig_no}_*600dpi.tif*",
                    f"Figure{fig_no}_*600dpi.tif*"
                ],
                f"Figure {fig_no} TIFF"
            )

    else:
        pdf = locate_unique_by_globs(
            FINAL_FIG_DIR,
            [
                f"Figure_{fig_no}_*.pdf",
                f"Figure{fig_no}_*.pdf"
            ],
            f"Figure {fig_no} PDF"
        )

        tif = locate_unique_by_globs(
            FINAL_FIG_DIR,
            [
                f"Figure_{fig_no}_*600dpi.tif*",
                f"Figure{fig_no}_*600dpi.tif*"
            ],
            f"Figure {fig_no} TIFF"
        )

    main_sources[fig_no] = {
        "pdf": pdf,
        "tiff": tif
    }

    print(f"Figure {fig_no}")
    print(f"  PDF  : {pdf}")
    print(f"  TIFF : {tif}")


# ============================================================
# 5. Supplementary figure paths
# ============================================================

SUPP_FIGURES = {

1: {
    "dir": SUPP_ROOT / "Supplementary_Figure_S1_gene_body_coverage",
    "pdf": "Supplementary_Figure_S1_gene_body_coverage.pdf",
    "tiff": "Supplementary_Figure_S1_gene_body_coverage_600dpi.tiff"
},

2: {
    "dir": SUPP_ROOT / "Supplementary_Figure_S2_coverage_sensitivity",
    "pdf": "Supplementary_Figure_S2_coverage_adjusted_sensitivity.pdf",
    "tiff": "Supplementary_Figure_S2_coverage_adjusted_sensitivity_600dpi.tiff"
},

3: {
    "dir": SUPP_ROOT / "Supplementary_Figure_S3_sequencing_mapping_assignment_QC",
    "pdf": "Supplementary_Figure_S3_sequencing_mapping_assignment_QC.pdf",
    "tiff": "Supplementary_Figure_S3_sequencing_mapping_assignment_QC_600dpi.tiff"
},

4: {
    "dir": SUPP_ROOT / "Supplementary_Figure_S4_fastp_QC_sensitivity",
    "pdf": "Supplementary_Figure_S4_fastp_QC_adjusted_expression_space_sensitivity.pdf",
    "tiff": "Supplementary_Figure_S4_fastp_QC_adjusted_expression_space_sensitivity_600dpi.tiff"
},

5: {
    "dir": SUPP_ROOT / "Supplementary_Figure_S5_2021_historical_reference_structure",
    "pdf": "Supplementary_Figure_S5_2021_historical_reference_expression_structure.pdf",
    "tiff": "Supplementary_Figure_S5_2021_historical_reference_expression_structure_600dpi.tiff"
},

6: {
    "dir": SUPP_ROOT / "Supplementary_Figure_S6_2021_DEG_overview",
    "pdf": "Supplementary_Figure_S6_2021_historical_reference_DEG_overview.pdf",
    "tiff": "Supplementary_Figure_S6_2021_historical_reference_DEG_overview_600dpi.tiff"
},

7: {
    "dir": SUPP_ROOT / "Supplementary_Figure_S7_cross_cohort_complete_diagnostics",
    "pdf": "Supplementary_Figure_S7_cross_cohort_complete_diagnostics.pdf",
    "tiff": "Supplementary_Figure_S7_cross_cohort_complete_diagnostics_600dpi.tiff"
},

8: {
    "dir": SUPP_ROOT / "Supplementary_Figure_S8_sample_exclusion_robustness",
    "pdf": "Supplementary_Figure_S8_sample_exclusion_and_leave_one_out_robustness.pdf",
    "tiff": "Supplementary_Figure_S8_sample_exclusion_and_leave_one_out_robustness_600dpi.tiff"
}
}


# ============================================================
# 6. Copy figures with canonical filenames
# ============================================================

banner("COPYING MAIN FIGURES")

for fig_no in range(1, 6):

    pdf_dst = MAIN_DIR / f"Figure_{fig_no}.pdf"
    tif_dst = MAIN_DIR / f"Figure_{fig_no}_600dpi.tiff"

    copy_file(
        main_sources[fig_no]["pdf"],
        pdf_dst
    )

    copy_file(
        main_sources[fig_no]["tiff"],
        tif_dst
    )

    print(
        f"Figure {fig_no}: PASS "
        f"[{image_info(tif_dst)}]"
    )


banner("COPYING SUPPLEMENTARY FIGURES")

for s in range(1, 9):

    info = SUPP_FIGURES[s]

    pdf_src = info["dir"] / info["pdf"]
    tif_src = info["dir"] / info["tiff"]

    if not pdf_src.exists():
        raise FileNotFoundError(
            f"Missing Supplementary Figure S{s} PDF:\n{pdf_src}"
        )

    if not tif_src.exists():
        raise FileNotFoundError(
            f"Missing Supplementary Figure S{s} TIFF:\n{tif_src}"
        )

    pdf_dst = (
        SUPP_FIG_DIR /
        f"Supplementary_Figure_S{s}.pdf"
    )

    tif_dst = (
        SUPP_FIG_DIR /
        f"Supplementary_Figure_S{s}_600dpi.tiff"
    )

    copy_file(
        pdf_src,
        pdf_dst
    )

    copy_file(
        tif_src,
        tif_dst
    )

    print(
        f"Supplementary Figure S{s}: PASS "
        f"[{image_info(tif_dst)}]"
    )


# ============================================================
# 7. Supplementary Tables
# ============================================================

banner("ASSEMBLING SUPPLEMENTARY TABLES")


# ------------------------------------------------------------
# Table S1: final 2026 metadata
# ------------------------------------------------------------

table_s1_src = recursive_find_unique(
    basenames=[
        "rattus_meta_2026_77samples.tsv",
        "02_2026_metadata_aligned.tsv"
    ],
    required=True
)

table_s1_dst = (
    SUPP_TABLE_DIR /
    "Supplementary_Table_S1_2026_primary_sample_metadata.tsv"
)

copy_file(
    table_s1_src,
    table_s1_dst
)

print("Supplementary Table S1: PASS")


# ------------------------------------------------------------
# Table S2: complete sequencing QC
# ------------------------------------------------------------

table_s2_src = (
    SUPP_ROOT /
    "Supplementary_Figure_S3_sequencing_mapping_assignment_QC" /
    "13B_fastp_QC" /
    "Supplementary_Table_S2_2026_primary_77sample_full_sequencing_QC.tsv"
)

if not table_s2_src.exists():
    table_s2_src = recursive_find_unique(
        basenames=[
            "Supplementary_Table_S2_2026_primary_77sample_full_sequencing_QC.tsv"
        ],
        required=True
    )

table_s2_dst = (
    SUPP_TABLE_DIR /
    "Supplementary_Table_S2_2026_primary_sequencing_QC.tsv"
)

copy_file(
    table_s2_src,
    table_s2_dst
)

print("Supplementary Table S2: PASS")


# ------------------------------------------------------------
# Table S3: 2026 timepoint-vs-Baseline DEG summary
# ------------------------------------------------------------

table_s3_rows = [
    {"Stage": "0 h",  "Up": 760,  "Down": 901,  "Total": 1661},
    {"Stage": "1 h",  "Up": 1344, "Down": 2389, "Total": 3733},
    {"Stage": "2 h",  "Up": 1088, "Down": 1237, "Total": 2325},
    {"Stage": "6 h",  "Up": 1158, "Down": 903,  "Total": 2061},
    {"Stage": "10 h", "Up": 477,  "Down": 977,  "Total": 1454},
    {"Stage": "14 h", "Up": 1169, "Down": 1435, "Total": 2604},
    {"Stage": "18 h", "Up": 1583, "Down": 1892, "Total": 3475},
    {"Stage": "36 h", "Up": 1207, "Down": 1422, "Total": 2629},
    {"Stage": "60 h", "Up": 1610, "Down": 1579, "Total": 3189},
    {"Stage": "72 h", "Up": 2230, "Down": 2224, "Total": 4454}
]

write_tsv(
    SUPP_TABLE_DIR /
    "Supplementary_Table_S3_2026_timepoint_vs_Baseline_DEG_summary.tsv",
    table_s3_rows,
    ["Stage", "Up", "Down", "Total"]
)

print("Supplementary Table S3: PASS")


# ------------------------------------------------------------
# Table S4: master temporal specificity annotation
# ------------------------------------------------------------

table_s4_src = recursive_find_unique(
    basenames=[
        "07_2026_master_temporal_specificity_annotation.tsv"
    ],
    regex=r"07_2026.*master.*temporal.*specificity.*annotation.*\.tsv$",
    required=True
)

table_s4_dst = (
    SUPP_TABLE_DIR /
    "Supplementary_Table_S4_2026_temporal_specificity_annotation.tsv"
)

copy_file(
    table_s4_src,
    table_s4_dst
)

print("Supplementary Table S4: PASS")


# ------------------------------------------------------------
# Table S5: final 2021 metadata
# ------------------------------------------------------------

table_s5_src = recursive_find_unique(
    basenames=[
        "02_2021_metadata_aligned.tsv",
        "rattus_meta_2021.tsv"
    ],
    required=True
)

table_s5_dst = (
    SUPP_TABLE_DIR /
    "Supplementary_Table_S5_2021_historical_reference_sample_metadata.tsv"
)

copy_file(
    table_s5_src,
    table_s5_dst
)

print("Supplementary Table S5: PASS")


# ------------------------------------------------------------
# Table S6: 2021 DEG summary
# ------------------------------------------------------------

table_s6_rows = [
    {"Stage": "4 h",  "Up": 2363, "Down": 1675, "Total": 4038},
    {"Stage": "8 h",  "Up": 2516, "Down": 2314, "Total": 4830},
    {"Stage": "12 h", "Up": 1729, "Down": 963,  "Total": 2692},
    {"Stage": "16 h", "Up": 1549, "Down": 1368, "Total": 2917},
    {"Stage": "20 h", "Up": 978,  "Down": 538,  "Total": 1516},
    {"Stage": "24 h", "Up": 1781, "Down": 1847, "Total": 3628},
    {"Stage": "48 h", "Up": 2173, "Down": 1403, "Total": 3576}
]

write_tsv(
    SUPP_TABLE_DIR /
    "Supplementary_Table_S6_2021_timepoint_vs_Baseline_DEG_summary.tsv",
    table_s6_rows,
    ["Stage", "Up", "Down", "Total"]
)

print("Supplementary Table S6: PASS")


# ------------------------------------------------------------
# Table S7: all 70 cross-cohort stage pairs
# ------------------------------------------------------------

table_s7_src = (
    SUPP_ROOT /
    "Supplementary_Figure_S7_cross_cohort_complete_diagnostics" /
    "Supplementary_Figure_S7_all_70_stage_pair_diagnostics.tsv"
)

if not table_s7_src.exists():
    table_s7_src = recursive_find_unique(
        basenames=[
            "Supplementary_Figure_S7_all_70_stage_pair_diagnostics.tsv"
        ],
        required=True
    )

table_s7_dst = (
    SUPP_TABLE_DIR /
    "Supplementary_Table_S7_cross_cohort_70_stage_pair_diagnostics.tsv"
)

copy_file(
    table_s7_src,
    table_s7_dst
)

print("Supplementary Table S7: PASS")


# ------------------------------------------------------------
# Table S8: combined sample-exclusion / leave-one-out summary
# ------------------------------------------------------------

s8_dir = (
    SUPP_ROOT /
    "Supplementary_Figure_S8_sample_exclusion_robustness"
)

s8_global_file = (
    s8_dir /
    "Supplementary_Figure_S8A_multi_sample_global_structure.tsv"
)

s8_cross_file = (
    s8_dir /
    "Supplementary_Figure_S8B_multi_sample_cross_cohort_robustness.tsv"
)

s8_loo_file = (
    s8_dir /
    "Supplementary_Figure_S8C_D_single_sample_leave_one_out.tsv"
)

for x in [
    s8_global_file,
    s8_cross_file,
    s8_loo_file
]:
    if not x.exists():
        raise FileNotFoundError(
            f"Missing S8 source table: {x}"
        )

global_rows = read_tsv(
    s8_global_file
)

cross_rows = read_tsv(
    s8_cross_file
)

loo_rows = read_tsv(
    s8_loo_file
)

cross_by_scenario = {
    r["Scenario"]: r
    for r in cross_rows
}

table_s8_rows = []

for r in global_rows:

    scenario = r["Scenario"]
    c = cross_by_scenario.get(
        scenario,
        {}
    )

    table_s8_rows.append({
        "Analysis_type": "Multi-sample exclusion",
        "Scenario_or_removed_sample": scenario,
        "N_samples": r.get("N_samples", ""),
        "PERMANOVA_R2": r.get("PERMANOVA_R2", ""),
        "PERMANOVA_P": r.get("PERMANOVA_P", ""),
        "PERMDISP_P": r.get("PERMDISP_P", ""),
        "Cross_cohort_matrix_correlation":
            c.get(
                "Correlation_with_primary_70cell_matrix",
                ""
            ),
        "Cross_cohort_matrix_MAE":
            c.get(
                "Mean_absolute_matrix_difference",
                ""
            ),
        "Minimum_genomewide_LFC_Spearman": "",
        "Review_type": ""
    })


for r in loo_rows:

    table_s8_rows.append({
        "Analysis_type": "Single-sample leave-one-out",
        "Scenario_or_removed_sample":
            r.get("Removed_sample", ""),
        "N_samples":
            r.get("N_samples", ""),
        "PERMANOVA_R2":
            r.get("PERMANOVA_R2", ""),
        "PERMANOVA_P":
            r.get("PERMANOVA_P", ""),
        "PERMDISP_P":
            r.get("PERMDISP_P", ""),
        "Cross_cohort_matrix_correlation":
            r.get(
                "Cross_cohort_matrix_correlation",
                ""
            ),
        "Cross_cohort_matrix_MAE":
            r.get(
                "Cross_cohort_matrix_MAE",
                ""
            ),
        "Minimum_genomewide_LFC_Spearman":
            r.get(
                "Minimum_genomewide_LFC_Spearman",
                ""
            ),
        "Review_type":
            r.get(
                "Review_type",
                ""
            )
    })


write_tsv(
    SUPP_TABLE_DIR /
    "Supplementary_Table_S8_sample_exclusion_and_leave_one_out_summary.tsv",
    table_s8_rows,
    [
        "Analysis_type",
        "Scenario_or_removed_sample",
        "N_samples",
        "PERMANOVA_R2",
        "PERMANOVA_P",
        "PERMDISP_P",
        "Cross_cohort_matrix_correlation",
        "Cross_cohort_matrix_MAE",
        "Minimum_genomewide_LFC_Spearman",
        "Review_type"
    ]
)

print("Supplementary Table S8: PASS")


# ============================================================
# 8. Supplementary Table titles / legends
# ============================================================

TABLE_LEGENDS = {

1: (
"Supplementary Table S1 | Sample metadata for the final 77-sample "
"2026 primary rat skeletal-muscle contusion cohort."
),

2: (
"Supplementary Table S2 | Per-sample sequencing, fastp quality, "
"HISAT2 alignment and featureCounts assignment metrics for the "
"final 77-sample 2026 primary cohort."
),

3: (
"Supplementary Table S3 | Numbers of up- and downregulated genes "
"for each 2026 post-injury stage relative to Baseline. DEGs were "
"defined using DESeq2 with adjusted P<0.05 and |log2 fold change|≥1."
),

4: (
"Supplementary Table S4 | Gene-level temporal annotations for the "
"2026 primary cohort, integrating continuous-time spline likelihood-"
"ratio-test results, amplitude filtering and TMM-normalized TPM-based "
"Tau, SPM and peak-stage annotations."
),

5: (
"Supplementary Table S5 | Sample metadata for the independently "
"reprocessed 60-sample 2021 historical-reference cohort."
),

6: (
"Supplementary Table S6 | Numbers of up- and downregulated genes "
"for each 2021 post-injury stage relative to its own Baseline group. "
"DEGs were defined using DESeq2 with adjusted P<0.05 and "
"|log2 fold change|≥1."
),

7: (
"Supplementary Table S7 | Complete diagnostics for all 70 cross-"
"cohort stage pairs, including Spearman and Pearson genome-wide "
"response concordance, shared-DEG counts and shared-DEG direction "
"concordance across 16,918 common independently filtered genes."
),

8: (
"Supplementary Table S8 | Multi-sample exclusion and single-sample "
"leave-one-out sensitivity results for global sampling-stage structure "
"and the 10 × 7 cross-cohort response-concordance matrix."
)
}


# ============================================================
# 9. Source Data collection
# ============================================================

banner("COLLECTING FIGURE SOURCE DATA")

source_index = []


def register_source(figure_id, src, description):
    src = Path(src)

    if not src.exists():
        return

    target_dir = SOURCE_DIR / figure_id
    target_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    target = (
        target_dir /
        safe_name(src.name)
    )

    copy_file(
        src,
        target
    )

    source_index.append({
        "Figure": figure_id,
        "Source_file": str(
            target.relative_to(PACKAGE)
        ),
        "Description": description,
        "Original_path": str(src)
    })


# ------------------------------------------------------------
# Main Figure 1 source
# ------------------------------------------------------------

register_source(
    "Figure_1",
    table_s1_src,
    "2026 final sample metadata used to define the primary cohort."
)

register_source(
    "Figure_1",
    table_s5_src,
    "2021 historical-reference sample metadata."
)


# ------------------------------------------------------------
# Main Figure 2-5 panel/source TSVs already generated
# ------------------------------------------------------------

panel_search_dirs = [
    FINAL_FIG_DIR,
    FINAL_FIG_DIR / "panels"
]

for fig_no in range(2, 6):

    found = []

    for d in panel_search_dirs:

        if not d.exists():
            continue

        found.extend(
            sorted(
                d.glob(
                    f"Figure_{fig_no}*.tsv"
                )
            )
        )

    for src in sorted(
        set(found)
    ):
        register_source(
            f"Figure_{fig_no}",
            src,
            f"Source data generated for main Figure {fig_no}."
        )


# ------------------------------------------------------------
# Supplementary Figure source TSVs
# ------------------------------------------------------------

for s in range(1, 9):

    d = SUPP_FIGURES[s]["dir"]

    for src in sorted(
        d.rglob("*.tsv")
    ):

        # Skip package-like tables already handled separately only
        # if desired; retaining them in Source Data is acceptable.
        register_source(
            f"Supplementary_Figure_S{s}",
            src,
            f"Source/diagnostic data for Supplementary Figure S{s}."
        )


# ============================================================
# 10. Write captions
# ============================================================

banner("WRITING MASTER CAPTIONS")


main_caption_file = (
    CAPTION_DIR /
    "Main_Figure_captions_Figure_1_to_5.txt"
)

with open(
    main_caption_file,
    "w"
) as f:

    for i in range(1, 6):
        f.write(
            MAIN_CAPTIONS[i].strip()
        )
        f.write(
            "\n\n"
        )


supp_caption_file = (
    CAPTION_DIR /
    "Supplementary_Figure_captions_S1_to_S8.txt"
)

with open(
    supp_caption_file,
    "w"
) as f:

    for i in range(1, 9):
        f.write(
            SUPP_CAPTIONS[i].strip()
        )
        f.write(
            "\n\n"
        )


all_caption_file = (
    CAPTION_DIR /
    "All_Figure_captions_master.txt"
)

with open(
    all_caption_file,
    "w"
) as f:

    f.write(
        "MAIN FIGURES\n"
    )
    f.write(
        "=" * 60 + "\n\n"
    )

    for i in range(1, 6):
        f.write(
            MAIN_CAPTIONS[i].strip()
        )
        f.write(
            "\n\n"
        )

    f.write(
        "\nSUPPLEMENTARY FIGURES\n"
    )
    f.write(
        "=" * 60 + "\n\n"
    )

    for i in range(1, 9):
        f.write(
            SUPP_CAPTIONS[i].strip()
        )
        f.write(
            "\n\n"
        )


table_legend_file = (
    CAPTION_DIR /
    "Supplementary_Table_titles_and_legends_S1_to_S8.txt"
)

with open(
    table_legend_file,
    "w"
) as f:

    for i in range(1, 9):
        f.write(
            TABLE_LEGENDS[i].strip()
        )
        f.write(
            "\n\n"
        )


print("Main Figure captions      : PASS")
print("Supplementary captions    : PASS")
print("Supplementary table legends: PASS")


# ============================================================
# 11. Caption index TSV
# ============================================================

caption_rows = []

for i in range(1, 6):
    caption_rows.append({
        "Type": "Main Figure",
        "Number": f"Figure {i}",
        "Caption": MAIN_CAPTIONS[i]
    })

for i in range(1, 9):
    caption_rows.append({
        "Type": "Supplementary Figure",
        "Number": f"Supplementary Figure S{i}",
        "Caption": SUPP_CAPTIONS[i]
    })

for i in range(1, 9):
    caption_rows.append({
        "Type": "Supplementary Table",
        "Number": f"Supplementary Table S{i}",
        "Caption": TABLE_LEGENDS[i]
    })


write_tsv(
    CAPTION_DIR /
    "Captions_and_legends_index.tsv",
    caption_rows,
    [
        "Type",
        "Number",
        "Caption"
    ]
)


# ============================================================
# 12. Source Data index
# ============================================================

write_tsv(
    CAPTION_DIR /
    "Source_Data_index.tsv",
    source_index,
    [
        "Figure",
        "Source_file",
        "Description",
        "Original_path"
    ]
)


# ============================================================
# 13. Supplementary Table index
# ============================================================

supp_table_rows = []

for i in range(1, 9):

    candidates = sorted(
        SUPP_TABLE_DIR.glob(
            f"Supplementary_Table_S{i}_*"
        )
    )

    if len(candidates) != 1:
        raise RuntimeError(
            f"Expected one Supplementary Table S{i}, "
            f"found {len(candidates)}."
        )

    supp_table_rows.append({
        "Table": f"Supplementary Table S{i}",
        "Filename": candidates[0].name,
        "Title_and_legend": TABLE_LEGENDS[i]
    })


write_tsv(
    CAPTION_DIR /
    "Supplementary_Table_index.tsv",
    supp_table_rows,
    [
        "Table",
        "Filename",
        "Title_and_legend"
    ]
)


# ============================================================
# 14. Figure file index
# ============================================================

figure_index = []

for i in range(1, 6):
    figure_index.append({
        "Figure": f"Figure {i}",
        "PDF": f"01_Main_Figures/Figure_{i}.pdf",
        "TIFF": f"01_Main_Figures/Figure_{i}_600dpi.tiff"
    })

for i in range(1, 9):
    figure_index.append({
        "Figure": f"Supplementary Figure S{i}",
        "PDF":
            f"02_Supplementary_Figures/"
            f"Supplementary_Figure_S{i}.pdf",
        "TIFF":
            f"02_Supplementary_Figures/"
            f"Supplementary_Figure_S{i}_600dpi.tiff"
    })


write_tsv(
    CAPTION_DIR /
    "Figure_file_index.tsv",
    figure_index,
    [
        "Figure",
        "PDF",
        "TIFF"
    ]
)


# ============================================================
# 15. README
# ============================================================

readme = """Scientific Data submission package

Directory structure
===================

01_Main_Figures
  Figure_1 to Figure_5.
  Each figure is provided as PDF and 600-dpi TIFF.

02_Supplementary_Figures
  Supplementary Figure S1 to Supplementary Figure S8.
  Each figure is provided as PDF and 600-dpi TIFF.

03_Supplementary_Tables
  Supplementary Table S1 to Supplementary Table S8.

04_Source_Data
  Figure-level source and diagnostic tables grouped by figure.

05_Captions_and_Indexes
  Master figure captions, supplementary-table legends,
  figure-file index and source-data index.

06_Manifests
  Complete package manifest and SHA256 checksums.

Analysis conventions
====================

2026 primary cohort:
  77 samples.
  Baseline = uninjured reference.
  R0h = immediate post-injury stage and is not Baseline.

2021 cohort:
  Previously published independently reprocessed historical-reference cohort.
  60 samples.

Cross-cohort analyses:
  Each cohort was analysed independently relative to its own Baseline.
  The cohorts were not pooled and no ComBat or cross-cohort batch
  correction was applied.

Sensitivity analyses:
  Adjusted expression matrices were used only as technical sensitivity
  analyses and do not replace the original primary expression matrices.

Review-flagged samples:
  Review flags were diagnostic and were not automatic exclusion criteria.
  All final 77 samples were retained in the primary 2026 analyses.

This package was assembled by:
  15A_assemble_submission_package.py
"""

with open(
    PACKAGE /
    "README_submission_package.txt",
    "w"
) as f:
    f.write(readme)


# ============================================================
# 16. Build package manifest
# ============================================================

banner("BUILDING MANIFEST / SHA256")

manifest_rows = []

for path in sorted(
    PACKAGE.rglob("*")
):

    if not path.is_file():
        continue

    if MANIFEST_DIR in path.parents:
        continue

    rel = path.relative_to(
        PACKAGE
    )

    manifest_rows.append({
        "Relative_path": str(rel),
        "Size_bytes": path.stat().st_size,
        "SHA256": sha256_file(path),
        "Image_metadata": image_info(path)
    })


manifest_file = (
    MANIFEST_DIR /
    "submission_package_manifest.tsv"
)

write_tsv(
    manifest_file,
    manifest_rows,
    [
        "Relative_path",
        "Size_bytes",
        "SHA256",
        "Image_metadata"
    ]
)


checksum_file = (
    MANIFEST_DIR /
    "SHA256SUMS.txt"
)

with open(
    checksum_file,
    "w"
) as f:

    for row in manifest_rows:
        f.write(
            f"{row['SHA256']}  "
            f"{row['Relative_path']}\n"
        )


# ============================================================
# 17. Final validation
# ============================================================

banner("FINAL VALIDATION")


# Main figures
for i in range(1, 6):

    pdf = (
        MAIN_DIR /
        f"Figure_{i}.pdf"
    )

    tif = (
        MAIN_DIR /
        f"Figure_{i}_600dpi.tiff"
    )

    if not pdf.exists() or not tif.exists():
        raise RuntimeError(
            f"Figure {i} incomplete."
        )


# Supplementary figures
for i in range(1, 9):

    pdf = (
        SUPP_FIG_DIR /
        f"Supplementary_Figure_S{i}.pdf"
    )

    tif = (
        SUPP_FIG_DIR /
        f"Supplementary_Figure_S{i}_600dpi.tiff"
    )

    if not pdf.exists() or not tif.exists():
        raise RuntimeError(
            f"Supplementary Figure S{i} incomplete."
        )


# Supplementary tables
for i in range(1, 9):

    hits = list(
        SUPP_TABLE_DIR.glob(
            f"Supplementary_Table_S{i}_*"
        )
    )

    if len(hits) != 1:
        raise RuntimeError(
            f"Supplementary Table S{i} validation failed."
        )


print("Main Figures 1-5             : PASS")
print("Supplementary Figures S1-S8 : PASS")
print("Supplementary Tables S1-S8  : PASS")
print(
    f"Source-data files            : {len(source_index)}"
)
print(
    f"Manifest files               : {len(manifest_rows)}"
)


# ============================================================
# 18. Final report
# ============================================================

banner("15A SUBMISSION PACKAGE COMPLETED")

print(f"Package:\n{PACKAGE}\n")

print("Primary caption files:")
print(
    "  05_Captions_and_Indexes/"
    "Main_Figure_captions_Figure_1_to_5.txt"
)
print(
    "  05_Captions_and_Indexes/"
    "Supplementary_Figure_captions_S1_to_S8.txt"
)
print(
    "  05_Captions_and_Indexes/"
    "Supplementary_Table_titles_and_legends_S1_to_S8.txt"
)

print()
print("Indexes:")
print(
    "  05_Captions_and_Indexes/"
    "Figure_file_index.tsv"
)
print(
    "  05_Captions_and_Indexes/"
    "Supplementary_Table_index.tsv"
)
print(
    "  05_Captions_and_Indexes/"
    "Source_Data_index.tsv"
)

print()
print("Manifest:")
print(
    "  06_Manifests/"
    "submission_package_manifest.tsv"
)
print(
    "  06_Manifests/"
    "SHA256SUMS.txt"
)

print()
print("PASS")
