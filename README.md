# mFLOW-Seq: Immunoglobulin class–resolved microbiota profiling across the female reproductive tract

[![Python](https://img.shields.io/badge/Python-3-blue.svg)](https://www.python.org/)
[![License: GPL v3](https://img.shields.io/badge/License-GPL_v3-blue.svg)](LICENSE)
[![Zenodo Version 6](https://img.shields.io/badge/Zenodo-Version_6-blue.svg)](https://doi.org/10.5281/zenodo.19437248)

Data, analysis scripts and figure-generation materials for the mFLOW-Seq study of immunoglobulin class–resolved microbiota across female reproductive tract compartments.

**Public archive: Version 6.** The submitted [GitHub repository](https://github.com/Prakashbio/mflow-seq-immune-microbiome-reproductive-tract) and [Zenodo concept DOI](https://doi.org/10.5281/zenodo.19437248) are unchanged. Start with [REPRODUCE.md](REPRODUCE.md) for the current analysis commands or the [Zenodo contents guide](https://zenodo.org/records/24505636/files/00_README_CONTENTS.md?download=1) for downloads. This version replaces the current publication workbook with the final 38-sheet version, which omits Table S2C.4. The historical participant-associated "Neg" reports are not established experimental negative controls and are not presented as such in the current publication workbook.

## Overview

This repository accompanies the manuscript:

> **Immunoglobulin class–resolved engagement of membrane-intact microbiota across female reproductive tract compartments using mFLOW-Seq**
>
> Prakash Lingasamy, Merli Saare, Sergio Vela Moreno, Karita Särekannu, Dmitri Lubenets, Alberto Sola-Leyva, Kasun Madhuranga Godakuma Godagedara, Naisarg Patel, Vijayachitra Modhukur, Reet Mändar, Andres Salumets
>
> *Manuscript under review.*

The retained sequencing dataset contains 149 libraries from 15 women and 40 recorded preparations. The current release includes the consolidated 38-sheet supplementary workbook, six main figures and nineteen numbered supplementary figures (23 supplementary graphics), their captions, and preserved analysis sources.

The README follows the original repository's section order. The public version is **6**. Internal build labels such as `v25.2.2` remain in original directory and script names to preserve working paths and provenance.

### About mFLOW-Seq

mFLOW-Seq combines fluorescence-activated cell sorting and 16S rRNA sequencing to profile PRE and immunoglobulin-associated microbial fractions from vaginal, cervical and endometrial preparations. The workflow links flow-cytometry records, sequencing-community profiles and participant/preparation-aware downstream analyses.

## Project Structure

```text
mflow-seq-immune-microbiome-reproductive-tract/
├── README.md                         # Project guide, retaining the original layout
├── REPRODUCE.md                      # Current reproduction instructions
├── verify_release_v25_2_2.py          # File integrity and source-link checks
├── replay_review_checks_v25_2_2.py    # Current numerical checks; optional figure replay
├── reproduce_public_figures_v25_2_2_RC2.py
├── SHA256SUMS.json                   # Current release payload hashes
├── REPRODUCIBILITY_SCOPE.json        # Reproduction entry-point information
├── LICENSE / LICENSE_SCOPE.txt / FILE_LICENSE_SCOPE.csv
├── revisions/
│   ├── revision_v25_2_2/
│   │   ├── selected_figures/         # Current figures and captions
│   │   ├── selected_tables/          # Current supplementary Excel workbook
│   │   └── Publication_table_source_index.csv
│   ├── revision_v25_2_1/             # Preserved previous revision outputs
│   ├── revision_v25/                 # Sequencing summaries, benchmarks and workflow files
│   └── revision_v24_8/
│       └── source_materials/        # Preserved analysis inputs, code and result records
├── CHANGELOG.md                     # Public release changes; earlier code stays in Git history
├── requirements.txt                  # Original core Python requirements
├── run_all.py / scripts/             # Retained historical analysis route
└── data/ / Output/                   # Retained historical inputs and outputs
```

The root replay scripts and `revisions/revision_v25_2_2/` are the current entry points. The earlier `run_all.py`, `scripts/`, `data/` and `Output/` files remain available with their original history.

## Features

- **Diversity:** Shannon, Simpson, observed richness and Pielou evenness.
- **Community composition:** taxon/genus summaries, Bray–Curtis ordination and matched sorted-to-own-PRE distances.
- **Repeated observations:** participant/preparation-aware analyses and saved resampling records.
- **Immunoglobulin fractions:** reconstructed FACS percentages and matched IgA/IgM sequencing representation summaries.
- **Storage and community groups:** descriptive comparisons using the recorded preparation metadata.
- **Sequencing records:** 149-library accounting, supplied benchmark tables and Nextflow/shell workflow files.
- **Verification:** file hashes, 80 current publication source blocks, numerical replay and 58 selected graphics/caption checks.

### Main Figures

Current files are in [selected_figures](revisions/revision_v25_2_2/selected_figures). The public figure driver rebuilds the selected set; each figure has a separate caption file.

| Current file | Figure | Description | Reproduction route |
|---|---|---|---|
| `Main_Figure_1.png` | 1 | mFLOW-Seq workflow overview | Public figure driver and preserved artwork |
| `Main_Figure_2.png` | 2 | PRE genus composition and observed woman/site profiles | Public figure driver |
| `Main_Figure_3.png` | 3 | FACS gate percentages and within-woman IgM-minus-IgA contrasts | Public figure driver |
| `Main_Figure_4.png` | 4 | Shannon diversity with preparation matching | Public figure driver |
| `Main_Figure_5.png` | 5 | Participant-connected PCoA and matched own-PRE distances | Public figure driver |
| `Main_Figure_6.png` | 6 | Matched IgA/IgM sequencing representation and relative abundances | Public figure driver |

---

### Supplementary Figures

| Current figure | Description | Current file(s) |
|---|---|---|
| S1 | Illustrative clinical flow-cytometry plots | `Supplementary_Figure_S1.emf` |
| S2 | Original FRT microscopy before/after sorting | `Supplementary_Figure_S2.emf` |
| S3 | Original saliva microscopy before/after sorting | `Supplementary_Figure_S3.emf` |
| S4 | Analytical feature-selection summaries | `Supplementary_Figure_S4.png` |
| S5 | PRE rank correlations by site | `Supplementary_Figure_S5.png` |
| S6 | Observed PRE assignment accumulation | `Supplementary_Figure_S6.png` |
| S7 | PRE profiles and retained assignment counts | `Supplementary_Figure_S7.png` |
| S8 | Across-site Shannon diversity within women | `Supplementary_Figure_S8.png` |
| S9 | Alpha-diversity metrics by site and fraction | `Supplementary_Figure_S9.png` |
| S10 | Matched sorted-to-own-PRE contrast intervals | `Supplementary_Figure_S10.png` |
| S11 | Paired Lactobacillus relative library abundances | `Supplementary_Figure_S11.png` |
| S12 | Genus-aggregated matched IgA/IgM representation | `Supplementary_Figure_S12.png` |
| S13 | PRE taxa by original vaginal-community description | `Supplementary_Figure_S13.png` |
| S14 | Alpha metrics by recorded storage | `Supplementary_Figure_S14.png` |
| S15 | Sequencing profiles by recorded storage | `Supplementary_Figure_S15.png` |
| S16 | PRE storage groups on common ordination axes | `Supplementary_Figure_S16.png` |
| S17 | Two endometrial preparations from participant 99 | `Supplementary_Figure_S17.png` |
| S18 | Retrospective HUT093 buffer-control reports | `Supplementary_Figure_S18_A/B/C.png` |
| S19 | Archived saliva controls and original pilot gates/comparator displays | `Supplementary_Figure_S19_A/B/C.png` |

S18 and S19 each contain three graphics; this gives 23 supplementary graphics across 19 numbered figures. Original microscopy artwork is unchanged. Full descriptions and panel definitions are in the accompanying caption files.

---

### Supplementary Data Scripts

| Script or source directory | Outputs | Description |
|---|---|---|
| `source_materials/own_PRE_between_sites_v24_5/` | Matched PRE distance summaries and verification results | Own-PRE numerical replay |
| `source_materials/score_provenance_bridge_v24_1/scripts/reproduce_score_bridge_v24_1.py` | Score bridge tables | Reconstructs 117 score links from the preserved sources |
| `source_materials/reviewer_metadata_v24_6/` | Participant/preparation/library metadata checks | Source metadata verification |
| `source_materials/baseline_v20_1_analysis_release/` | Preserved numerical inputs, models and resampling records | Baseline source bundle and Q15 checks |
| `revisions/revision_v25/scripts/recompute_v25_sequencing_summaries.py` | Sequencing-accounting/benchmark summaries | Recomputes arithmetic from supplied source tables |

Paths beginning `source_materials/` are relative to `revisions/revision_v24_8/`. The current [80-block source index](revisions/revision_v25_2_2/Publication_table_source_index.csv) links publication tables to preserved source files and hashes. The two source-index entries for withdrawn Table S2C.4 are excluded; other block identifiers are retained without renumbering.

---

### Shared Utilities (`shared_utils.py`)

The original `scripts/shared_utils.py` is retained for the historical scripts. It contains helpers for loading the original count/metadata tables, diversity metrics, transformations, distances and ordination. Current reproduction is coordinated by the root replay driver, which records a log and exit status for each step.

| Current utility | Purpose |
|---|---|
| `verify_release_v25_2_2.py` | Check payload hashes, source links and declared file membership |
| `replay_review_checks_v25_2_2.py` | Run the existing numerical checks and optional figure replay |
| `reproduce_public_figures_v25_2_2_RC2.py` | Rebuild the 58 selected graphics/captions using preserved inputs and artwork |

## Requirements

### System Requirements

- Python 3 for the verifier and numerical/figure drivers.
- R and `Rscript` for the figure route; the recorded plotting environment uses R 4.4.1.
- Poppler's `pdftoppm` on PATH for source-plot rendering.
- Storage for the extracted source files and a separate output directory.

### Python Dependencies

The original core requirements are retained in `requirements.txt`:

```text
numpy >= 1.23
pandas >= 1.5
matplotlib >= 3.6
scipy >= 1.10
scikit-learn >= 1.2
```

Source-specific environment records are retained under `revisions/revision_v24_8/source_materials/baseline_v20_1_analysis_release/environment/`. The plotting session record is in `revisions/revision_v24_8/source_materials/current_v24_figure_reproduction/audit/R_sessionInfo.txt`; it identifies ggplot2, patchwork, vegan and their recorded dependencies. The saliva-pilot source directory also retains its original requirements file. Consult these records and the selected source script when running individual historical analyses.

## Installation

### 1. Clone or Download the Repository

```bash
git clone https://github.com/Prakashbio/mflow-seq-immune-microbiome-reproductive-tract
cd mflow-seq-immune-microbiome-reproductive-tract
```

Alternatively, use [Zenodo Version 6](https://doi.org/10.5281/zenodo.24505636). Download either `mFLOW_seq_Version6_complete.zip` or all six numbered topic ZIPs. Extract every topic ZIP into the same empty directory, preserving stored paths. The topic archives partition the same files as the complete ZIP; the exact payload count is given in the contents guide and file index. Keep the standalone download guide and index outside that extracted analysis directory.

### 2. Create a Virtual Environment

Create the environment outside the checkout so the integrity verifier sees only the declared source files.

**Windows (PowerShell):**

```powershell
python -m venv ../.mflow-seq-venv
..\.mflow-seq-venv\Scripts\Activate.ps1
```

**Linux/macOS:**

```bash
python3 -m venv ../.mflow-seq-venv
source ../.mflow-seq-venv/bin/activate
```

### 3. Install Dependencies

```bash
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

For the public numerical and figure replay, `requirements-replay-tested.txt` records the Python package versions used in the verified replay, including Pillow, pypdf, reportlab, openpyxl and lxml in addition to the original core packages. Install the source-specific R dependencies and Poppler required by the selected route, using the environment records above. [REPRODUCE.md](REPRODUCE.md) gives the current commands.

## Usage

### Run All Analyses

Run integrity and existing numerical checks from the extracted root:

```bash
python3 verify_release_v25_2_2.py
python3 replay_review_checks_v25_2_2.py --output /absolute/path/to/new_checks
```

To include the existing figure replay:

```bash
python3 replay_review_checks_v25_2_2.py --output /absolute/path/to/new_checks_with_figures --figures
```

Use a different, previously nonexistent output directory outside the repository for each run. The replay writes step logs and `Replay_receipt.json`.

### Run Specific Figures

The current driver reproduces the selected figures and captions together:

```bash
python3 reproduce_public_figures_v25_2_2_RC2.py --output /absolute/path/to/new_figure_output
```

It writes `figures/` and `Public_figure_replay_receipt.json` under the chosen output directory. The earlier `run_all.py --fig ...` command remains associated with the historical figure numbering.

### List Available Scripts

Current script paths and publication outputs are indexed in the tables above and in the source index. The original script listing remains available for the historical analysis route:

```bash
python3 run_all.py --list
```

### Run Individual Scripts

Example: replay the supplied sequencing-summary arithmetic into a new external folder:

```bash
python3 revisions/revision_v25/scripts/recompute_v25_sequencing_summaries.py --input revisions/revision_v25/data --output /absolute/path/to/new_sequencing_summary
```

Individual source scripts retain their own parameters and original environment records.

### Command-Line Options

```text
python3 replay_review_checks_v25_2_2.py --output NEW_EXTERNAL_DIRECTORY [--figures]

--output    Required: a new directory outside the extracted repository
--figures   Also run the current figure-reproduction route
-h, --help  Show command help
```

The standalone figure driver also requires `--output NEW_EXTERNAL_DIRECTORY`. Existing results are not overwritten by these drivers.

## Output Files

Current publication outputs are preserved under `revisions/revision_v25_2_2/`. New replay outputs are written to the external directory selected with `--output`.

### Main Figures

`selected_figures/Main_Figure_1.png` through `Main_Figure_6.png`, each with its corresponding `_caption.txt` file.

### Supplementary Figures

`selected_figures/Supplementary_Figure_S1.emf` through S3 retain the original artwork. S4–S17 use PNG files; S18 and S19 each use A/B/C PNG panels. All have caption files. The selected set totals 29 graphics plus 29 captions.

### Supplementary Data Tables and Statistical Outputs

- `selected_tables/Supplementary_Tables_CLEAN.xlsx`: final consolidated 38-sheet workbook; Table S2C.4 is not included.
- `Publication_table_source_index.csv`: 80 current publication source blocks with source paths and hashes.
- `Replay_receipt.json`: numerical replay status and per-step logs.
- `Public_figure_replay_receipt.json`: selected graphics/caption reproduction status.
- Source CSV/TSV, RDS/NPZ and environment records remain beside their original analysis directories.

# Data Files

### Input Data Description

| File or location | Description |
|---|---|
| `revisions/revision_v25_2_2/selected_tables/Supplementary_Tables_CLEAN.xlsx` | Current consolidated publication workbook |
| `revisions/revision_v25_2_2/Publication_table_source_index.csv` | Publication-table to source-file crosswalk |
| `revisions/revision_v25/data/` | Supplied library accounting, sequencing QC and computational benchmark tables |
| `revisions/revision_v24_8/source_materials/` | Coded count/abundance matrices, metadata, recorded FACS summaries and preserved analytical sources |
| `SHA256SUMS.json` | SHA-256 for every listed release payload |
| `EXISTING_GITHUB_FILES_RETAINED.json` | Hashes of 53 retained historical GitHub files outside the standalone revision ZIP |
| `data/otu_table.csv`, `sample_metadata.csv`, `taxonomy.csv` | Retained original public input tables |
| `data/facs_binding_events.csv` and other original `data/` files | Historical FACS/event and analysis inputs retained unchanged |

The public release includes exactly three sample-free HUT093 buffer FCS files. The 34 human-sample pilot FCS files remain private. Public pilot-display reproduction uses the included saved aggregate bins and plotting parameters.

**Data availability:** [Submitted Zenodo concept DOI](https://doi.org/10.5281/zenodo.19437248); [specific public Version 6](https://doi.org/10.5281/zenodo.24505636); [contents guide](https://zenodo.org/records/24505636/files/00_README_CONTENTS.md?download=1). The complete ZIP is supplied alongside six topic ZIPs, a file index and download checksums. Replacing the current publication workbook does not erase prior Git commits or published Zenodo versions; those historical snapshots can still contain withdrawn Table S2C.4.

Raw sequencing records are available under NCBI BioProject PRJNA1454428 / study SRP692419. The supplied workflow records identify near-full-length 27F/1492R amplification, Tn5-fragmented libraries, NextSeq 1000 single-end reads delivered at 80 nt, Kraken2 Standard and Bracken `-r 75`; the latter is the Bracken reference-distribution setting.

## Analysis Methods

### Statistical Tests

- Shannon-diversity models accounting for repeated participant/preparation observations, with descriptive Simpson, richness and evenness summaries.
- Bray–Curtis dissimilarities calculated from library relative abundances, with recorded restricted-permutation PERMANOVA/dispersion analyses.
- Matched own-PRE comparisons, including between-site contrasts and saved resampling records.
- Site-specific IgA/IgM representation summaries with participant-aware aggregation and saved bootstrap intervals.
- Recorded storage/community-group summaries and supplied sequencing-QC/benchmark arithmetic.

Exact model outputs, settings, resampling objects and source scripts are retained with their analysis records. Figure captions and the source index identify the current outputs.

### Color Scheme

Use the keys embedded in the current figures: colour and symbol assignments follow each figure's site, fraction or participant legend. The original scientific figure colours and microscopy are retained. Publication tables use black text and white backgrounds.

## Troubleshooting

### Common Issues

**1. Missing Python/R package or external command**

Activate the intended Python environment, install the core requirements, and check the source-specific Python/R environment records. Ensure `Rscript` and `pdftoppm` are on PATH for the figure route. Inspect the failed step's log in the external output directory.

**2. File not found**

Extract the complete ZIP, or all six topic ZIPs, into one empty directory. Preserve stored paths and run the root driver from that directory. Use the current entry points rather than the older `run_all.py` for current outputs.

**3. Existing or internal output directory**

Choose a new output directory outside the extracted repository. The current drivers reject existing output directories and do not overwrite prior results.

**4. Zenodo ZIP preview message**

The complete archive contains more entries than Zenodo's ZIP browser preview displays. Download/extract it directly, or use the six smaller topic ZIPs. The readable contents guide is the default preview. Both extraction routes reconstruct the same verified files.

**5. Unlisted source file during integrity checks**

Keep replay outputs and the standalone download guide/index outside the extracted analysis folder. The verifier accepts the declared release files, plus the 53 recorded legacy GitHub files when present.

## Contact

For questions, issues or collaboration inquiries:

- **Prakash Lingasamy:** prakash.lingasamy@ut.ee
- **Andres Salumets:** andres.salumets@ki.se

## License

Project analysis code retains the GNU General Public License v3.0 in [LICENSE](LICENSE). Project-owned coded tabular data and derived scientific outputs use the dataset record's CC BY 4.0 terms. Source-specific notices and attribution remain applicable; [LICENSE_SCOPE.txt](LICENSE_SCOPE.txt) and [FILE_LICENSE_SCOPE.csv](FILE_LICENSE_SCOPE.csv) describe the file categories.

## Citation

If you use these materials, cite the manuscript and the public data/code archive:

> Lingasamy P, Saare M, Vela Moreno S, Särekannu K, Lubenets D, Sola-Leyva A, Godagedara KMG, Patel N, Modhukur V, Mändar R, Salumets A.
> **Immunoglobulin class–resolved engagement of membrane-intact microbiota across female reproductive tract compartments using mFLOW-Seq**
> *Manuscript under review.*

Archive citation: Lingasamy, P., & Godakumara, K. (2026). *Immunoglobulin class–resolved engagement of membrane-intact microbiota across female reproductive tract compartments using mFLOW-Seq* (Version 6) [Dataset]. Zenodo. https://doi.org/10.5281/zenodo.24505636

The submitted concept DOI, https://doi.org/10.5281/zenodo.19437248, continues to resolve to the latest version.

Version 6 replaces the selected supplementary workbook with the final submission version and refreshes the current source index, release checks and download manifests. The final workbook includes the later table/report clarifications and corrections, not merely a sheet deletion. Preserved analytical inputs, scripts and figure artwork are not being reclassified or claimed to have been newly experimentally validated. Earlier Zenodo versions remain historical snapshots rather than the current submission workbook.

## Acknowledgments

We thank the study contributors and research participants. The original repository acknowledgment of support from the Estonian Research Council and the Horizon Europe / European Research Council framework is retained.
