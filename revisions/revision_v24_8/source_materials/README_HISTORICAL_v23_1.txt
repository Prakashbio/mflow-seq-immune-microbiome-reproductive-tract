mFLOW-Seq v23.1 prepared GitHub/Zenodo update — 11 November 2026

LOCAL STAGING ONLY: nothing has been uploaded or published.

metadata_corrections contains the current HUT093 buffer-control annotations.
Read CONTROL_METADATA_AMENDMENT.txt first. It supersedes historical statements
that the named HUT093 BSA/PBS1/PBS2 tubes contained antibodies or DAPI.

current_v23_figure_reproduction contains the unchanged current panel figures,
captions, vector exports, plotting inputs and scripts. Reproduce into a new
folder using its existing instructions:
  python3 current_v23_figure_reproduction/scripts/reproduce_v23_figures.py current_v23_figure_reproduction NEW_OUTPUT_FOLDER

baseline_v20_1_analysis_release preserves all 580 files from the original
public-update baseline, including analysis inputs, scripts and saved records.
Historical code and data are retained unchanged. The current figure folder
supersedes the baseline's figure organization; metadata_corrections supersedes
earlier control-composition annotations. All prior scientific payload bytes are
unchanged. No new models, tests or re-gating were performed in v23.1.

The figure driver was independently run for v23 and reproduced all 23 selected
PNG/EMF assets exactly. This downstream reproducibility does not validate the
upstream sequencing protocol, species accuracy, DAPI selection or antibody
specificity. Buffer-only controls do not supply antibody-reagent background
measurements. Do not describe this staging folder as a published release.

Private author correspondence, AI conversations, reviewer reports and manuscript
drafts have not been added to this metadata amendment. The historical guide is
retained as README_HISTORICAL_v23.txt for provenance.
