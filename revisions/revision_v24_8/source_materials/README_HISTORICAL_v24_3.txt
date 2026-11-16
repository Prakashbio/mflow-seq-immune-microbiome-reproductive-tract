mFLOW-Seq v24.3 prepared GitHub/Zenodo update — 12 November 2026

LOCAL STAGING. No upload or publication has been performed.

V24.2 incorporates the laboratory-confirmed workflow, proprietary-method limits,
and independent-extraction/shared-batch mock design. It retains the v24.1 score
bridge and adds the existing-record 149-library metadata crosswalk. Counts, metadata, primary/sensitivity results, figures and microscopy
are unchanged from v24. Historical versions of the README and manifests are
retained with HISTORICAL in their names. All other existing files retain their
exact bytes.

REPRODUCE THE SCORE BRIDGE
python3 score_provenance_bridge_v24_1/scripts/reproduce_score_bridge_v24_1.py THIS_RELEASE_ROOT NEW_OUTPUT_FOLDER

The driver verifies seven packaged input-file hashes and all four output hashes.
It independently reconstructs all 149 source-vector correspondences and 117
primary site/taxon scores while separating source-factor allocation, matched
pair eligibility and woman weighting. This is a provenance diagnostic, not a
new statistical model or new primary estimate. Original sample codes are kept.
Its source-vector trace, original script and verification record are supplied.
The upstream FASTQ identity and execution limitations remain unverified.

REPRODUCE THE CURRENT FIGURES (unchanged from v24)
python3 current_v24_figure_reproduction/scripts/reproduce_v24_figures.py THIS_RELEASE_ROOT NEW_OUTPUT_FOLDER

The current driver retains the v23.3 observed-PRE Figure 2 layout and horizontal
Figure 5A keys, with the v24 label/legend refinements to 17 plots. All 29 current
image assets are hash checked. Original microscopy/FACS artwork, colours,
symbols, groups, intervals and plotted values are preserved. The v24.3 document
changes do not alter figure captions or image files, so no additional plotting
version or model refit is needed. The previous independent 29-asset build and
layer-preservation records are retained.

The full numerical baseline and downstream scripts remain under
baseline_v20_1_analysis_release. The source-coded plotting inputs remain under
current_v23_figure_reproduction/source_data. Earlier software/environment and
reproduction records retain their original versions. These materials reproduce
the documented downstream analyses; they do not authenticate missing laboratory
protocol history, taxonomic accuracy, raw-read processing or DAPI validation.

The current revised Word/PDF documents, private coauthor message and reviewer
reports are supplied in the separate author package, not this repository tree.

Laboratory-confirmed primer/amplification/read descriptions do not authenticate original raw-read processing. Exact proprietary formulation limits are disclosed in the current Methods. Metadata in metadata_update_149_2026_11_12 has explicit source-stage and missing-data distinctions. No public SRA correction or repository upload was performed.

The facs_control_retrospective_2026_11_12 directory adds three original no-antibody buffer PDFs, transparent population transcriptions, historical-versus-retrospective comparisons, unchanged report renders and a regeneration script. These are retrospective analyses of the same historical acquisitions; no new experiment, clinical correction or statistical reanalysis. The original source workbooks, statistical code/data and score bridge remain unchanged.
