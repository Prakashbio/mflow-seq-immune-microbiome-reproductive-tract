mFLOW-Seq v23.2 prepared GitHub/Zenodo update — 11 November 2026

LOCAL STAGING ONLY. Nothing has been uploaded or published.

CURRENT FIGURE ORGANIZATION AND REPRODUCTION
current_v23_2_figure_reproduction contains current captions, expected graphic
hashes, numbering crosswalks and the current driver. Six main and 19 supplementary
figures are selected; S18 and S19 each have three image plates (29 total assets).
Rebuild into a new folder from this release root:
  python3 current_v23_2_figure_reproduction/scripts/reproduce_v23_2_figures.py THIS_RELEASE_ROOT NEW_OUTPUT_FOLDER

The driver rebuilds the unchanged 23 v23 graphics from their frozen inputs,
renumbers the supplementary files and regenerates six added control/pilot plates.
No analytical model is refitted. Requirements: the recorded R environment for the
old figure driver, Python numpy/Pillow/pypdf/reportlab and Poppler pdftoppm.

facs_control_pilot_display_v23_2 contains 37 unmodified source FCS files, original
operator graph EMFs, original HUT093 report, display scripts, captions and source
hash/crop/bin records. Read its README for the distinction between original gates
and retrospective, ungated exported-event displays.

CURRENT METADATA
metadata_corrections/Table_S2V_5B_FACS_control_types_UPDATED_v23_2_2026-11-11.csv
is the current control-type annotation. Earlier metadata files/guides remain
historical. The HUT093 buffers contain no antibodies or DAPI. Early unstained and
IgA/IgM-only saliva records are available, without establishing clinical coverage.

UNCHANGED ANALYTICAL BASELINE
baseline_v20_1_analysis_release retains all original numerical-analysis inputs,
scripts, saved draws and tests. current_v23_figure_reproduction retains the prior
figure build used by the new driver. Do not interpret old file numbering or
historical guides as the current submission selection.

VALIDATION LIMITS
Graph reproduction does not authenticate upstream FASTQ processing, resolve the
reported-versus-registered sequencing protocol, validate species assignments,
establish quantitative DAPI concordance or demonstrate antibody specificity.
No private author correspondence, reviewer report or manuscript draft is added
to this new public staging content. No external release action was performed.
