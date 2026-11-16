mFLOW-Seq v24 prepared GitHub and Zenodo update — 11 November 2026

LOCAL STAGING ONLY. No upload or publication has been performed.

CURRENT FIGURE REPRODUCTION
python3 current_v24_figure_reproduction/scripts/reproduce_v24_figures.py THIS_RELEASE_ROOT NEW_OUTPUT_FOLDER

The driver reconstructs the preserved v23.3 baseline, including the observed-PRE
Figure 2 layout and Figure 5A horizontal legends, then applies the v24 layout
script to 17 plots. All 29 selected image assets are checked against their current
hashes. Source values, points, intervals, box summaries, heatmap values, colours,
fraction symbols, groups, scale transformations and statistical results are
retained. S15 annotations omit the redundant n= prefix, with every displayed
count separately checked identical. L. abbreviates Lactobacillus; Fr/Fz denote
fresh/frozen on selected axes. Full taxon/source identifiers remain unchanged.

current_v24_figure_reproduction holds the 17 revised PNG/PDF/SVG exports, scripts,
current captions, data-layer signatures, count-label audit and figure crosswalk.
The 12 other assets retain their exact bytes, including original microscopy,
FACS composite/control/pilot graphics, workflow and current Figures 2 and 5.

The unchanged frozen source inputs are under current_v23_figure_reproduction/
source_data and baseline_v20_1_analysis_release. The layout script copies inputs
to temporary storage before evaluating the historical plotting code. No source
result is overwritten and no new statistical model is fitted. Earlier scripts
and guides remain historical records. Use the current v24 driver for the current
file selection and captions. The recorded R and Python/Poppler dependencies are
retained in the earlier guides and session records; no new library is required.

Six undefined one-taxon-library Pielou values remain missing as before; no values
are filled or discarded for this layout update. Downstream reproducibility does
not establish the unresolved sequencing/read-processing provenance, taxonomic
accuracy, quantitative DAPI agreement, clinical antibody background or viability.
No manuscript, author correspondence or reviewer text is added to this staging.
