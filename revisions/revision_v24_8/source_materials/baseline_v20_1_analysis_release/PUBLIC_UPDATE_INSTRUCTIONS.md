# Prepared v20.1 repository and deposition update

This folder is the complete local downstream update payload. Nothing has been
pushed, uploaded or submitted. Keep the previous repository/deposit history.
Place this payload in a versioned revision directory when making the later
authorized repository update; do not overwrite historical root analyses blindly.

Use README.md and REPRODUCE_v20_1.md for supported replay commands. The existing
reproduced/ files are frozen reference outputs included in SHA256SUMS.json and
must be included in the commit/deposit. The model refit creates additional files;
run it in a separate working copy. Figure and joint runners require a new output
directory outside this candidate. Transient Python/R files are ignored.

The original GPL v3 code license is retained. No new data license, release tag,
commit identifier, DOI, publication date or ORCID has been invented. The planning
metadata identifies the manuscript creators and existing record lineage. Resolve
actual release/deposition fields at the time of the authorized update. This is
not an API upload payload or a completed publication claim.

This payload excludes the manuscript/rebuttal drafts, private correspondence,
recruitment records, raw FCS acquisitions and raw FASTQ files. Its successful
downstream replays do not reproduce the upstream FASTQ processing or establish
the missing direct quantitative DAPI validation.

The consolidation includes filter reconciliation, a read-only Q15 verifier and documented provenance aliases. No numerical input or statistical result has changed.
