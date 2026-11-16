Existing FACS control and pilot evidence for v23.2

These are original experimental files and source graphs, not new validation
experiments. HUT093 BSA/PBS1/PBS2 contain buffer only without antibodies or DAPI.
The original report regions are software gates. No original gates were moved.

Figures S18A–C preserve cropped native DAPI/SSC and IgA/IgM source graph panels.
Figure S19A displays all 12 early saliva unstained/single-colour acquisitions
from exported events; S19B shows representative original DAPI pilot gates;
S19C displays six candidate SYTO9/PI records. The retrospective bin displays use
all recorded events, arcsinh(value/150), 150 bins per axis and log colour counts;
they apply no additional compensation or reconstructed sorting/live/dead gate.
Acquisition counts do not establish independent biological n or paired aliquots.

Reproduce into a new folder:
  python3 scripts/plot_v23_2_facs_evidence.py THIS_FOLDER NEW_OUTPUT_FOLDER

Requires Python with numpy, Pillow, pypdf and reportlab, base R/Rscript and
Poppler pdftoppm on PATH. No matplotlib, X11 or cairo installation is required.
The native graph bitmaps are extracted from the supplied operator EMFs; source
PDF/EMF/FCS files are checked unchanged before/after plotting. The source report
contains all original scatter, IgG and population-table panels as well.

This does not establish quantitative DAPI/comparator agreement, original
compensation history, matched clinical FMO/background coverage, per-cell binding
specificity, viability or post-sort purity. All controls, including PBS2, remain.
No background subtraction or new inferential test is performed.
