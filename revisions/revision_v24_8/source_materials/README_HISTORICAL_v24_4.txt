mFLOW-Seq v24.4 prepared repository update — 12 November 2026

Existing participant-aware analyses, source-coded data, recorded permutations,
bootstrap draws, figure scripts and score-provenance bridge are retained.
No clinical data, model or original microscopy changed in this update.
This archive does not reproduce the unresolved upstream FASTQ-to-table chain.

Current buffer figure: S18A–C is in current_v24_4_control_figures/. It uses the
operator-supplied BSA/PBS1/PBS2 PDFs exported 12 November 2026 from the same
historical acquisitions. Source report pixels/gates were not edited. Earlier
plots remain in historical_reports_NOT_CURRENT_FIGURES/. Both original and
retrospective counts are documented in the paper supplement. There is no
clinical background subtraction or antibody-specific false-positive estimate.
The S20 filenames retained in facs_control_retrospective_2026_11_12 are labels
from historical working v24.3; the current paper uses these reports as S18.

To verify the control arithmetic and original PDF hashes:
python3 facs_control_retrospective_2026_11_12/verify_control_statistics.py
To regenerate current control report renders (requires Poppler pdftoppm):
python3 current_v24_4_control_figures/reproduce_current_S18.py
To replay all 117 primary site/taxon scores without changing inputs:
python3 score_provenance_bridge_v24_1/scripts/reproduce_score_bridge_v24_1.py . NEW_OUTPUT_DIRECTORY

Baseline release folders and README_HISTORICAL files preserve earlier scripts
and descriptions. Current numerical-analysis definitions are documented in the
participant-aware baseline release and current figure reproduction records;
historical equivalence/causal claims are not current conclusions.

Laboratory-confirmed near-full-length PCR/85-bp read layout does not authenticate
original run records, filtering, database/classifier settings or usable raw
access. Proprietary formulation limits and the restricted validation scope
are disclosed in the manuscript/rebuttal. Public metadata corrections and
GitHub/Zenodo upload have not been performed. This is local staging only.
