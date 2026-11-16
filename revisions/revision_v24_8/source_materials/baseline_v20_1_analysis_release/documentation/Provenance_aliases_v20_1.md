# Portable provenance aliases

PROJECT_ROOT denotes the original mFLOW-Seq project root. EXPERIMENT_ROOT denotes
the external source-experiment root. WORKSTATION_ROOT denotes other historical
user-root references. These are provenance aliases, not promised bundled files.
Original input basenames, worksheet/cell locators and source hashes are retained.
The prior sealed v20 archive preserves the full original provenance strings.

Two archived Python preparation scripts now accept MFLOW_NATIVE_WORKBOOK for
their external source workbook; the default is that workbook basename relative
to the working directory. Their other historical project-layout conventions
remain. Use the supported replay runners for recorded downstream reproduction;
these two historical scripts are not claimed as complete standalone pipelines.
No data values, participant codes, phase or storage fields were removed.
