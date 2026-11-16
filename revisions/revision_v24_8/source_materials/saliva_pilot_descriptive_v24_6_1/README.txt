Descriptive historical saliva pilot evidence

Original source FCS files remain unchanged in ../facs_control_pilot_display_v23_2/source_fcs. The 13 June 2023 table compares all-event exported detector medians in numbered Ig-labelled and unstained records. The 20 June table records candidate same-number acquisitions and their start-time gaps. Indices are not authenticated donors or independently paired biological samples. No gate, additional compensation, p-value, confidence interval, specificity or DAPI/SYTO9–PI agreement was inferred. Detector voltage, gain, range/amplification and trigger/SPILL fields matched within the reported contrasts. Source hashes and complete settings accompany the tables.

Medians and counts reproduce the supplied independent summaries. Small percentile differences reproduce float32 probability/interpolation arithmetic; the independent table uses float64 quantiles. This does not change the reported median contrasts.

Run scripts/verify_pilot_summaries.py with --source pointing to the original source_fcs folder and --out to an empty results folder; --provided-folder can optionally point to the three supplied comparison CSVs.

Reproduction from this folder:
python scripts/verify_pilot_summaries.py --source ../facs_control_pilot_display_v23_2/source_fcs --out replay

The recorded numerical environment is NumPy 2.3.5 (requirements.txt). No random number generator is used. Original FCS files are read only.
