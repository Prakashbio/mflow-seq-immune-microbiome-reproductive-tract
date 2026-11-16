# mFLOW-Seq reviewer revision v24.8

Current revised downstream analyses, source-coded data and figure reproduction for **Immunoglobulin class-resolved microbial profiling across female reproductive tract compartments using mFLOW-Seq**.

Start in [revisions/revision_v24_8](revisions/revision_v24_8/REPRODUCE.md). The versioned directory contains the source materials, current figure/caption selection, one consolidated supplementary workbook and read-only verification scripts. It represents the scientific data/figures used in the audited v24.8 revision; v24.8 standardises publication table presentation, removes internal prose, corrects retrospective control captions and clarifies the detector-panel rationale. Scientific source data and all figure pixels remain those of the audited v24.6.1 staging.

The retained dataset contains 149 libraries from 15 women and 40 preparations. Original source identifiers/counts and microscopy are retained. Analyses account for participant/preparation dependence, include the joint site × fraction sensitivity and within/across-site own-PRE contrasts, and report site-specific relative-representation scores and appropriate exploratory limitations.

The original root `run_all.py`, `scripts/`, `data/` and `Output/` belong to the historical public version. They are preserved in repository history and are not the entry point for current revised figures. The former README is retained in `historical_documentation/README_before_revision_at_5f85149.md`; its historical claims do not describe the current inference.

## Evidence scope and access

DAPI-negative selection is operational and was selected after early SYTO9/PI panel interference. No new dye-comparison experiment or quantitative agreement result is claimed. Buffer-only controls are instrument/background observations, not antibody-binding specificity measurements. Species labels remain provisional. The data do not establish biological ascension or preservation equivalence.

The laboratory reports 27F/1492R V1–V9 amplification, tagmentation and nominal 85-bp single-end sequencing. Original FASTQ execution/QC provenance is not supplied by this downstream archive. As checked 12 November 2026, SRA BioProject PRJNA1454428 contains 149 unreleased run records scheduled for 1 June 2027; its reviewer metadata link does not provide FASTQ downloads. Registered MiSeq/V3–V4 descriptions are inconsistent with the laboratory clarification. Read access and metadata corrections require separately authorized action. No private reviewer token is included here.

This is a prepared local update; no new repository commit/tag, Zenodo DOI or public availability is asserted by this file. Historical Zenodo concept DOI: https://doi.org/10.5281/zenodo.19437248 (observed version v2: https://doi.org/10.5281/zenodo.19449798).

## License

The existing GPLv3 code license is retained in LICENSE. Existing data/source licenses and attribution are preserved; this preparation assigns no new data license. Original source artwork and microscopy are retained with their provenance.

## Confirmed HUT093 control composition

The author confirmed on 13 November 2026 that the three named BSA, PBS1 and PBS2 acquisitions contained DAPI but no immunoglobulin-detection antibodies or biological sample. The selected workbook, S18 captions and amended control-type CSV reflect this reagent metadata correction. All historical and retrospective counts, clinical results, statistical calculations and figure pixels remain unchanged. The 1466 original source files are preserved byte-for-byte as historical records; their superseded reagent annotations must not be treated as the current control description. The source index directs D066 to the versioned amendment. This is prepared locally, not an online release.
