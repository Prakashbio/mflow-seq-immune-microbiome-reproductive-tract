# mFLOW-Seq revision analysis preparation

This local candidate collects the source-coded inputs, revised analytical outputs, saved randomization records, working figures and scientific scripts associated with the Microbiome major revision. It is synchronized to the consolidated v20.1 author-review selection prepared 6 November 2026 and has not been uploaded. The recorded experimental and upstream-processing limitations remain explicit. No upload has occurred.

The experimenter confirms 17 distinct recruits and 16 women contributing to the sequencing workflow. The available clinical/FACS records concern 16 women; these stages are kept distinct. The retained sequencing dataset comprises 149 libraries from 15 women and 40 preparations. Original library codes are preserved. Species labels remain provisional; the quantitative DAPI-concordance request and parts of the upstream sequencing workflow remain unresolved.

Start with `REPRODUCE_v20_1.md` for the current read-only checks and replay commands. Earlier dated records remain historical evidence.

## Files

- `bundle/` retains source-coded count matrices, metadata, revised tables, saved bootstrap/permutation/model objects and the current working figure assets/captions.
- `analysis_scripts/source_archive/` retains the recorded scientific scripts with their historical layout conventions and stable root aliases; two workbook paths accept MFLOW_NATIVE_WORKBOOK. These are provenance sources, not a claim of portable execution.
- `analysis_scripts/refit_recorded_models.R` changes only three root/output path assignments in the previously verified fresh refit. It checks 10 Shannon, four focused beta and three dispersion models by replaying recorded permutations.
- `verify_bundle.py` checks the manifest and independently recalculates the primary site-score arithmetic and selected FACS/mock arithmetic using Python's standard library.
- `documentation/` provides file hashes, source provenance, table columns, current figure mapping and the direct submitted/current comparison.
- `environment/` retains the recorded R/Python environments. These are not a complete installation lock for every archived figure script.

## Verification

From this directory:

```sh
python3 verify_bundle.py
Rscript --vanilla analysis_scripts/refit_recorded_models.R .
```

The R command requires the packages and compatible versions listed in `environment/Fresh_R_refit_sessionInfo.txt`. It writes new files under `reproduced/fresh_model_check/`. Manifest/source verification does not recreate the raw FASTQ pipeline or validate species calls or DAPI. The portable figure runner below staged the recorded scripts and inputs in an isolated layout: eight plotting scripts reproduced all 19 statistical PNG figures with identical bytes and pixels in the recorded local environment. A fresh-machine installation and the upstream FASTQ pipeline remain unverified.

## Interpretation and provenance

The primary IgA/IgM scores are ratios of equal-woman means within site after averaging repeated preparations, with pointwise bootstrap intervals. They do not measure per-cell coating or affinity. Fraction/site R² estimates from different conditional models must not be ranked. Storage comparisons are descriptive and do not establish equivalence. The original extraction-blank and mock summaries have bounded interpretation. Original microscopy image bytes are preserved.

The 23 registered current figures and their synchronized captions are included; three original native EMF supplementary panels have now been displayed and checked in a local native Word export; Prakash manages final author details and EndNote references. Recheck all labels after any subsequent manuscript edit. The original repository's `run_all.py` and figure numbering describe historical analyses and must not be presented as regenerating these revised results. The repository history should retain the original version while the revision receives its own documented release.

## License and citation

The existing repository code-license text is preserved in `LICENSE`. Final data licensing and full release citation metadata must be reconciled with the existing Zenodo record and final author list; no new license or DOI is assigned in this preparation. Cite the finalized version-specific DOI and commit/tag once created.

## v16 additions

Table S2F.7 inventories 16 labelled SYTO9/PI pilot acquisitions. Table S2C.4 and the accompanying CSVs/original reports describe 15 historical filename-selected candidate controls. They are not established as 15 matched extraction blanks. The participant-99 storage interpretation is updated to confirmed same-biopsy fragments; numerical model outputs are unchanged. Redundant numbered copies are omitted where byte identity or the explicit current-figure registry establishes the canonical selection.

This remains a local candidate. Historical scripts and source provenance retain original workstation paths documented in `documentation/Local_path_review.csv`. The supplied model, figure and joint runners provide the tested execution paths; retaining a historical script does not mean it is a portable entry point. No private correspondence or author-review Word document is selected.

## Targeted v17 changes and figure replay

Numerical inputs, saved draws, statistical results and image bytes are unchanged from v16. The S2 caption now reflects the experimenter-confirmed clinical preparation linkage and largely dried, coverslipped microscopy. Sequencing location and the historical Mageeibacillus exclusion rationale are documented in the Word revision. Those text clarifications are not new experimental validation.

To regenerate the 19 statistical images from this candidate using the recorded R dependencies and Python Pillow, choose a new output directory outside the candidate:

```sh
python3 analysis_scripts/replay_recorded_figures.py --candidate . --output ../figure_replay --rscript Rscript
```

The four nonstatistical assets are preserved rather than regenerated. The completed local replay records are in `documentation/figure_reproduction/`. This runner does not reprocess FASTQ reads, establish DAPI accuracy or make every historical script portable.

## v18 response-location correction

The scientific inputs, numerical outputs, scripts and figure bytes are unchanged from v17. The reviewer response now cites native Word page/line locations from eight clean documents. This is a document cross-reference correction, not a new analysis or experimental validation.

## v19 documentation update

The numerical analytical payload, saved randomization records, scripts and figures are unchanged from v18. The manuscript now retains the author-confirmed recruitment total and explicitly acknowledges historical no-sample IgM-region occupancies as a specificity limitation. Main/diversity native Word locations were refreshed; the 29 responses include 229 supporting passage links and 138 grouped references. These document checks do not provide quantitative DAPI validation or reproduce FASTQ processing. The supplementary crosswalk and logical component inventory are supplied at the package root.

## v20 joint factorial addition

The revised manuscript, supplementary information and all 29 reviewer responses
now include the additional joint site × fraction inferential sensitivity.
Tables S2R.8–S2R.10 report its primary terms, three related sensitivity arms
and adjacent dispersion checks. The existing 149-library descriptive partition
and focused conditional tests remain distinct. All 23 figure images are
unchanged; one phrase in the S12 caption was clarified consistently.

`bundle/joint_factorial_q15/` contains original source-coded input matrices,
the 124-library eligibility manifest, 120-profile crosswalk, relative abundances,
results, model matrices, all 9,999 permutation rows and null statistics,
1,000 woman-cluster bootstrap draws per arm, model objects, session information
and independent numerical verification. Both HUT099 endometrial preparations
contribute through equal averaging within fraction. Original sample labels
are preserved. The source cohort still contains 149 libraries from 15 women;
the added complete-fraction model concerns 14 women and 30 anatomical units.

Run the full joint analysis and independent verification into a NEW directory
outside this candidate:

```sh
python3 analysis_scripts/run_joint_factorial_q15.py --candidate . --output ../joint_q15_reproduced --rscript Rscript --verify
```

R 4.4.1 with vegan 2.7-2 and permute 0.9-8 was used; full session information
is recorded in `environment/Joint_factorial_Q15_R_sessionInfo.txt`.
The joint runner and verifier require Python 3 with NumPy and pandas.
The custom Gower projection engine holds nuisance terms fixed and implements
Freedman–Lane reduced-model residual permutations explicitly. It does not
substitute vegan's partial-model permutation routine for this fixed-design
scheme. Primary term p-values use Holm across three tests; dispersion uses a
separate two-test family. Related sensitivities are not independent replication.
Fraction dispersion differs. Inference is approximate under the specified
residual-exchangeability assumptions; numerical agreement does not establish
those assumptions, isolated centroid effects or antibody causality.

The current public entry points are the root README and REPRODUCE.md. Earlier numerical source files and the original code licence are retained.

## v20.1 consolidation on 6 November 2026

The short supplementary filter clarification and `bundle/filter_reconciliation/`
record the six historical No/retained cases. All six pass the current original-code
rule; no count matrix, feature membership, model result or image changed.
The Q15 verifier now reads its analysis directory without writing into it by
default. Use --report-out with a new external JSON path if a saved report is wanted.
Historical workstation-root strings have been replaced by documented provenance
aliases; two archived external workbook paths accept MFLOW_NATIVE_WORKBOOK.
This cleanup does not claim that all historical scripts are standalone or that a
fresh-machine installation or upstream FASTQ workflow has been reproduced.
The scientific data-sharing/license scope still follows the authorized deposit.
The original source inputs and numerical baseline remain preserved.
