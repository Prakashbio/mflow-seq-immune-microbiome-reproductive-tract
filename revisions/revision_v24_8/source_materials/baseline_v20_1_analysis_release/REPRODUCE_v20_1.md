# Reproduce the consolidated v20.1 update

Run these commands from this candidate directory. The first three are read-only.

1. Integrity and selected downstream arithmetic:

   `python3 verify_bundle.py`

2. Reconstruct all 41 tracking flags, 287 original global cells and the six
   disputed retained cases using source-code and legacy-heading groups:

   `python3 analysis_scripts/verify_filter_reconciliation.py`

3. Independently verify the stored Q15 calculation without changing the input:

   `python3 analysis_scripts/verify_Q15_joint_factorial_2026_11_05.py bundle/joint_factorial_q15`

   Optionally pass `--report-out ../NEW_q15_verification.json`. The report path
   must be new and outside the analysis directory. Input-tree targets and
   existing reports are rejected. NumPy and pandas are required; SciPy is not.

4. Regenerate the 19 statistical figures into a new external output directory:

   `python3 analysis_scripts/replay_recorded_figures.py --candidate . --output ../NEW_figure_replay --rscript Rscript`

5. Refit the 17 recorded models IN A SEPARATE WORKING COPY, since this command
   writes new files below reproduced/fresh_model_check:

   `Rscript --vanilla analysis_scripts/refit_recorded_models.R .`

6. Recompute all joint-model sensitivity arms with 9,999 permutations and 1,000
   cluster bootstrap draws per arm into a NEW external directory:

   `python3 analysis_scripts/run_joint_factorial_q15.py --candidate . --output ../NEW_joint_replay --rscript Rscript --verify`

   The optional verifier writes its report next to that new directory, outside
   the analysis tree. Existing outputs/reports are preserved.

Use the recorded environments for R dependencies. R 4.4.1, vegan 2.7-2 and
permute 0.9-8 were used for the joint model. Figure replay additionally needs
Python Pillow and the R packages listed in the saved sessions. This documents
local recorded-environment replay, not a clean-machine installation test.

Historical source scripts retain historical project-layout conventions. Two
external workbook paths now accept MFLOW_NATIVE_WORKBOOK; see the provenance
alias documentation. They are not the supported portable entry points above.
The original sealed v20 source remains available for provenance comparison.

This payload does not supply raw FASTQs/FCS, authenticate the upstream classifier
execution, establish quantitative DAPI concordance or resolve data-sharing terms.
Source-coded participant, phase and storage metadata are preserved for analysis.
No repository push, Zenodo upload or new version DOI is asserted.
