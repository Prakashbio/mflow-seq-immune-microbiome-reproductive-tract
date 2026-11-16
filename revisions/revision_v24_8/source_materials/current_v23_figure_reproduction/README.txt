Current integrated v23 figures — 11 November 2026

From the parent folder, run:
  python3 current_v23_figure_reproduction/scripts/reproduce_v23_figures.py current_v23_figure_reproduction NEW_OUTPUT_FOLDER

Dependencies: Rscript with ggplot2, patchwork and svglite; Python with numpy, pandas, Pillow, pypdf and lxml. Package versions are recorded in audit/R_sessionInfo.txt.

Inputs are unchanged stored CSV/RDS results, not new model fits. The wrapper reproduces all 23 current PNG/EMF assets exactly, retaining microscopy and mapping the old standalone own-PRE plot into Main Figure 5B. S13 now shows the saved four alpha-diversity indices. The historical driver initially renders both intermediate displays before the wrapper applies the current mapping. No new statistical tests are run.

The supporting CSV Optional_S18_plotted_values.csv is the historical source filename for the current S13 data. Its 596 cells (149 libraries × four metrics) include six undefined Pielou values. It is not an additional current Figure S18.

Reproducibility here means the documented downstream figures and saved-result calculations. It does not establish bacterial viability, antibody-binding specificity, an authenticated upstream FASTQ workflow or validated species resolution. Repository upload has not been performed.
