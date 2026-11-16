# Reproduce and verify revision v25.2.2

Run from the extracted repository root. Python3 alone is sufficient for integrity and supplied-summary arithmetic. Numerical replay additionally requires the recorded NumPy/pandas/scipy environment. Figure replay also requires Rscript with the recorded R packages and Poppler pdftoppm on PATH. No fresh-machine lock or complete proprietary FASTQ-to-table replay is claimed.

## Integrity

    python3 verify_release_v25_2_2.py

This verifies the root manifest, all82 indexed source hashes and all1633 members retained from the preceding code archive. Run it before and after any replay.

## Existing numerical checks

    python3 replay_review_checks_v25_2_2.py --output /absolute/path/to/NEW_mflow_checks

The output must be new and outside the repository. This tested wrapper executes baseline/hash checks, the saved Q15 numerical check, reviewer metadata check,117-score bridge, own-PRE replay plus its independent verifier, and the1000 sequencing-summary arithmetic checks. It saves per-step logs and a receipt. It does not overwrite source data or fit a new clinical model.

## Optional full figure replay

    python3 replay_review_checks_v25_2_2.py --output /absolute/path/to/NEW_mflow_figures --figures

This also runs the existing figure-reproduction wrapper and checks all58 selected graphics/captions against stored hashes. Microscopy and native source artwork are preserved. Graphics reproducibility is environment-dependent; use the recorded R/Python sessions. Successful replay does not validate antibody specificity, DAPI agreement, upstream sequencing provenance or observational exchangeability.

## Standalone own-PRE independent verifier

The correct script path is:

    revisions/revision_v24_8/source_materials/own_PRE_between_sites_v24_5/verify_own_PRE_between_sites_v24_5.py

Use --workdir on a separate directory containing inputs/ and results/. The replay wrapper prepares that structure safely. The shortened command in the historical revision_v24_8 guide omitted this source path; it is superseded by this current guide.

## Current files

Use revisions/revision_v25_2_2/selected_tables/Supplementary_Tables_CLEAN.xlsx for the39-sheet publication workbook and revisions/revision_v25_2_2/selected_figures/ for the current figures/captions. Publication_table_source_index.csv links82 source blocks to unchanged versioned CSVs. Sequencing source summaries/main.nf/run.sh remain in revisions/revision_v25. Older source directories preserve their own historical guides and selected files; do not use their earlier workbook as the current supplement.

The supplied Nextflow fragments require external execution configuration, samples, databases and proprietary processing not fully supplied here. Simulation-generator/accession/seed records and other experimental source limits remain disclosed. No online access or new version DOI is asserted by these commands.
