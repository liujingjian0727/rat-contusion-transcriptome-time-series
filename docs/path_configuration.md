# Local path configuration

Some technical-validation scripts require access to sequencing-QC files
that are not distributed through this GitHub repository.

The following scripts use the environment variable `SCIDATA_DATA_ROOT`:

- `supplementary/12A_make_Supplementary_Figure_S1_gene_body_coverage.R`
- `supplementary/13A_make_Supplementary_Figure_S3_sequencing_mapping_assignment_QC.R`
- `supplementary/13B_integrate_fastp_QC_and_stage_diagnostic.R`

`SCIDATA_DATA_ROOT` should point to the local directory corresponding to
`bulk_RNA_2026`.

For example:

    export SCIDATA_DATA_ROOT=/path/to/bulk_RNA_2026

The expected local structure is:

    bulk_RNA_2026/
    ├── GeneBodyCoverage/
    ├── CleanData/
    │   └── fastq_stats.txt
    ├── hisat2/
    ├── run_featurecounts/
    └── fastp/
        └── fastp_QC_Cohort_2026.tsv

The public scripts do not require the original server-specific absolute path.

Core count, TPM and metadata inputs are described separately in
`docs/input_files.md`.
