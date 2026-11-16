# Reproduce the current downstream revision

Run from `revisions/revision_v24_8`. Use a new external output directory for each replay. Read-only verification requires Python 3; the numerical modules additionally need NumPy/pandas and the recorded R environments. Figure replay uses Python NumPy/pandas/Pillow/pypdf/lxml, Rscript with the stored packages and Poppler pdftoppm. Recorded sessions are retained in source_materials. This is not a fresh-machine installation lock or FASTQ-to-table pipeline.

1. Verify every current release file: `python3 verify_release.py`.
2. Check baseline data and selected arithmetic: `python3 source_materials/baseline_v20_1_analysis_release/verify_bundle.py`.
3. Check the saved joint model: `python3 source_materials/baseline_v20_1_analysis_release/analysis_scripts/verify_Q15_joint_factorial_2026_11_05.py source_materials/baseline_v20_1_analysis_release/bundle/joint_factorial_q15`.
4. Check the reviewer-facing metadata: `python3 source_materials/reviewer_metadata_v24_6/verify_reviewer_metadata_v24_6.py source_materials/reviewer_metadata_v24_6`.
5. Rebuild the 117-score bridge: `python3 source_materials/score_provenance_bridge_v24_1/scripts/reproduce_score_bridge_v24_1.py source_materials ../../NEW_score_bridge`.
6. Rebuild the cross-site own-PRE results: `python3 source_materials/own_PRE_between_sites_v24_5/own_PRE_between_sites_v24_5.py --inputs source_materials/own_PRE_between_sites_v24_5/inputs --output ../../NEW_own_PRE`. Its independent verifier writes a report, so run it only in a separate working copy: `python3 verify_own_PRE_between_sites_v24_5.py --workdir WORKING_COPY`.
7. Rebuild current figures/captions: `python3 reproduce_current_figures.py --output ../../NEW_current_figures`. The wrapper preserves source imagery, renders the updated BSA/PBS reports and checks the final selected asset hashes, including the current Figure 5 caption.
8. Reinspect the historical saliva detector summaries: `python3 source_materials/saliva_pilot_descriptive_v24_6_1/scripts/verify_pilot_summaries.py --source source_materials/facs_control_pilot_display_v23_2/source_fcs --out ../../NEW_saliva_summaries`.

The baseline `REPRODUCE_v20_1.md` describes optional refits and full joint-model reruns. Run these in a separate working copy, because some historical verification/refit helpers write new receipts. Historical working paths retained in provenance scripts are not supported portable entry points.

Current figures/tables are selected_figures/ and selected_tables/. Intermediate historical figure names and earlier captions in source_materials are provenance, not the journal selection. No raw-read access, contaminant status, species accuracy, viability or matched clinical antibody-specificity result is inferred from a successful replay.

## Publication table index

`Publication_table_source_index.csv` maps all 78 neutral workbook block identifiers D001–D078 to existing source CSVs and SHA256 digests. Paths are relative to this version directory, not assertions of an online release. The selected workbook contains the unchanged measurements with black-and-white publication styling and neutral prose; earlier workbooks in source_materials remain historical records. No participant or library code was renamed.

## Confirmed HUT093 control composition

The author confirmed on 13 November 2026 that the three named BSA, PBS1 and PBS2 acquisitions contained DAPI but no immunoglobulin-detection antibodies or biological sample. The selected workbook, S18 captions and amended control-type CSV reflect this reagent metadata correction. All historical and retrospective counts, clinical results, statistical calculations and figure pixels remain unchanged. The 1466 original source files are preserved byte-for-byte as historical records; their superseded reagent annotations must not be treated as the current control description. The source index directs D066 to the versioned amendment. This is prepared locally, not an online release.
