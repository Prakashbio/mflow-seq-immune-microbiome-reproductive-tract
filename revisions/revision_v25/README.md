# Sequencing QC and benchmark source files

This addition preserves the previous downstream scripts, input matrices, saved resampling records and selected figures exactly. The original archive payloads are retained under their existing paths. New material is under `revisions/revision_v25`.

## Added material
- Original supplied main.nf and example run.sh (not executed here).
- Original per-library QC workbook and lossless 149-row CSV export.
- Original single-species, equal-mixture and pair-mixture numerical CSVs.
- A portable Python standard-library arithmetic check, its outputs and code correspondence tables.

Run after extraction:

    python3 revisions/revision_v25/scripts/recompute_v25_sequencing_summaries.py --input revisions/revision_v25/data --output v25_summary_check_output

The script verifies supplied-summary balances and recalculates denominators, mixture complements and rank preservation. It does not generate simulations, run a classifier or rerun clinical analyses. The direct-versus-clade definition differences are flagged, not corrected. The source abs_error_pp column is signed; a name does not change its numerical definition.

Original workflow paths require inputs.json, samples, database paths and container parameters not provided here. The original proprietary application/list matching and historical filtering-stage bridge are incomplete; set_inputs.py and simulation-generation sources were mentioned but not supplied. A fully replayable upstream workflow is not claimed.

The new QC values belong to their labelled stages. Their totals do not equal archived report roots, and one 78EP/78Epre2 correspondence is unconfirmed. The clinical matrices, participant structure, original microscopy and existing figure/data analysis remain unchanged. The three simulation designs are computational benchmarks, separate from the physical D6310 mock extractions. Species labels remain provisional.

These files are included in the published supporting archive. Use the root README and REPRODUCE.md for the current public entry points. The supplied run.sh is a historical launcher example: its workflow/input paths refer to the original execution environment, not a complete runnable distribution.
