# Reproduce and verify revision v25.2.1

Run these commands from the repository root after extraction. Use a new external output directory for each replay. Python 3 alone suffices for the root hash check and sequencing-summary arithmetic. Scientific replay additionally needs the Python/R packages documented in the preserved session/environment records. This is not a fresh-machine installation lock or a complete upstream FASTQ pipeline.

## Check the release

    python3 verify_release_v25_2_1.py

This verifies every file against SHA256SUMS.json, all82 indexed source hashes, all1560 preserved source-archive members and selected figure/workbook hashes. It does not modify scientific inputs or infer experimental validation.

## Recompute the supplied sequencing-summary arithmetic

    python3 revisions/revision_v25/scripts/recompute_v25_sequencing_summaries.py --input revisions/revision_v25/data --output NEW_v25_sequencing_check

This checks1000 supplied-summary identities and recalculates denominators/complements/ranks. It does not run Kraken2/Bracken or generate simulations.

## Downstream numerical and figure replay

The exact existing commands for the baseline bundle, saved Q15 model, 117-score bridge, all nine between-site own-PRE comparisons and current figure reproduction are retained in:

    revisions/revision_v24_8/REPRODUCE.md

Follow its stated working directory, and write replay outputs to a new external directory. Historical selected_tables in revision_v24_8 predate Sergio's additions. Use revision_v25_2_1/selected_tables for the39-sheet workbook and revision_v25_2_1/selected_figures for the current manuscript selection. Figure pixels/captions match the preceding release; the current source index additionally links the four new source tables.

Original software, reference/database configuration, proprietary processing and simulation-generation gaps remain described in the manuscript and the versioned README files. A successful arithmetic replay does not resolve those gaps or validate observational exchangeability.
