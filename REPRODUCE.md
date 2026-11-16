# Public revision v25.2.2 — release candidate2

Run from the extracted public archive root. Requires the recorded Python/R packages and Poppler for figure rendering. No fresh-machine dependency lock or complete upstream pipeline is claimed.

    python3 verify_release_v25_2_2.py
    python3 replay_review_checks_v25_2_2.py --output /absolute/path/to/NEW_checks
    python3 replay_review_checks_v25_2_2.py --output /absolute/path/to/NEW_figures --figures

Use a new output directory outside the repository each time. Existing downstream numerical checks and analysis inputs are unchanged. The optional figure step uses `reproduce_public_figures_v25_2_2_RC2.py`; it patches a temporary external source copy to route the existing figure chain to the new public-summary renderer. Canonical source scripts are retained unchanged. All 58 selected graphics/captions must match saved hashes.

34 human-sample pilot FCS files are private and absent from this public release. S19A/C displays are re-rendered from the exact saved aggregated bins and verified rendering parameters. This reproduces displays; it does not reanalyse private event arrays,verify donor pairing,validate DAPI or establish clinical antibody specificity. Historical raw-pilot entry points require private inputs and are not the public reproduction route.

Current39-sheet workbook/figures are in `revisions/revision_v25_2_2/`;the 82-block source index links retained underlying tables. Sequencing QC/benchmarks are in `revisions/revision_v25/`. The supplied Nextflow fragments still require external execution configuration,reference databases and proprietary processing not fully supplied. Simulation generation/accession/seed information remains incomplete. Retain the78EP/78Epre2 correspondence caveat;do not relabel input libraries.

The archive is a revision overlay,not a full clone of the existing GitHub repository. Existing repository files remain online;the public root verifier covers exactly this archive’s listed files and explicit private-file exclusions.

In a full GitHub checkout,the root verifier also accepts and validates the 53 unchanged legacy files listed in EXISTING_GITHUB_FILES_RETAINED.json. Their Git blob hashes must match;other unlisted source files remain errors. Those legacy files are intentionally absent from the standalone revision ZIP.
