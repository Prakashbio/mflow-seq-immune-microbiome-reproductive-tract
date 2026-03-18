# mFLOW-Seq Microbiome Analysis Pipeline

[![Python](https://img.shields.io/badge/Python-3.9+-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-Research-orange.svg)](LICENSE)

Bioinformatics analysis pipeline and figure generation scripts for the mFLOW-Seq microbiome study investigating spatial immunoglobulin engagement of live microbiota across female reproductive tract compartments.

## Overview

This repository contains the complete computational analysis workflow used to generate all main and supplementary figures for the manuscript:

> **Spatial Immunoglobulin Engagement of Live Microbiota Across Female Reproductive Tract Compartments Revealed by mFLOW-Seq**
>
> Prakash Lingasamy, Merli Saare, Kasun Madhuranga Godakuma Godagedara, Sergio Vela Moreno, Karita Särekannu, Dmitri Lubenets, Alberto Sola-Leyva, Naisarg Patel, Vijayachitra Modhukur, Reet Mändar, Andres Salumets
>
> *Manuscript under communication*

### About mFLOW-Seq

mFLOW-Seq (Microbiome Flow-Sorting and Sequencing) is a novel methodology combining fluorescence-activated cell sorting (FACS) and 16S rRNA sequencing to profile immunoglobulin-bound bacterial populations in complex microbiome samples. This approach enables the characterization of host-microbiota interactions at the single-cell level.

## Project Structure

```
mFLOW-Seq/
│
├── run_all.py              # Master pipeline script
├── requirements.txt        # Python dependencies
├── README.md              # This file
│
├── data/                  # Input data files
│   ├── otu_table.csv             # OTU abundance matrix
│   ├── sample_metadata.csv       # Sample experimental metadata
│   ├── taxonomy.csv              # Taxonomic assignments
│   ├── facs_binding_events.csv   # FACS binding data
│   └── ...                       # Additional analysis tables
│
├── scripts/               # Analysis and figure generation scripts
│   ├── shared_utils.py           # Common functions and utilities
│   ├── fig2_alpha_diversity.py   # Main figure 2
│   ├── fig3_beta_diversity.py    # Main figure 3
│   ├── ...                       # Additional figure scripts
│   └── sfig*.py                  # Supplementary figure scripts
│
└── new_figures/           # Output directory (generated)
    ├── Figure2.png
    ├── Figure3.png
    └── ...
```

## Features

The pipeline performs comprehensive microbiome analyses including:

- **Alpha Diversity**: Shannon, Simpson indices, and species richness
- **Beta Diversity**: PCoA with Bray-Curtis dissimilarity and PERMANOVA
- **Differential Abundance**: Statistical testing for compositional differences
- **Taxonomic Composition**: Phylum and genus-level profiling
- **Immunoglobulin Binding**: FACS-based characterization (IgA, IgG, IgM)
- **Sample Comparison**: Fresh vs. frozen, anatomical sites, health conditions
- **Network Analysis**: Co-occurrence and prevalence networks
- **Quality Control**: Contamination analysis and read depth assessment

## Requirements

### System Requirements

- **Operating System**: Windows 10/11, Linux, or macOS
- **Python**: Version 3.9 or higher
- **RAM**: Minimum 8 GB recommended
- **Storage**: ~500 MB for data and outputs

### Python Dependencies

```
numpy >= 1.23
pandas >= 1.5
matplotlib >= 3.6
scipy >= 1.10
scikit-learn >= 1.2
```

## Installation

### 1. Clone or Download the Repository

```bash
git clone <repository-url>
cd final_scripts_1803
```

Or download and extract the ZIP archive.

### 2. Create a Virtual Environment

**Windows (PowerShell)**:
```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
```

**Linux/macOS**:
```bash
python3 -m venv .venv
source .venv/bin/activate
```

### 3. Install Dependencies

```bash
pip install --upgrade pip
pip install -r requirements.txt
```

## Usage

### Run All Analyses

Execute the complete pipeline to generate all figures:

```bash
python run_all.py
```

Expected runtime: ~5-15 minutes depending on system specifications.

### Run Specific Figures

Generate only selected main figures (e.g., Figures 2, 4, 8):

```bash
python run_all.py --fig 2 4 8
```

Generate specific supplementary figures:

```bash
python run_all.py --fig S1 S2 S10
```

### List Available Scripts

View all figure scripts and their outputs:

```bash
python run_all.py --list
```

### Run Individual Scripts

Execute a single script directly from the repository root:

```bash
python scripts/fig2_alpha_diversity.py
```

**Note**: Always run scripts from the repository root directory to maintain correct relative paths.

### Command-Line Options

```
python run_all.py [OPTIONS]

Options:
  --fig NUM [NUM ...]   Generate specific figures (e.g., 2 4 S1 S5)
  --list                List all available figure scripts
  -q, --quiet           Suppress detailed execution output
  -h, --help            Show help message
```

## Output Files

All generated files are saved to the `new_figures/` directory:

### Main Figures
- `Figure2.png` - Alpha diversity analysis
- `Figure3.png` - Beta diversity and PCoA
- `Figure4.png` - FACS binding events
- `Figure5.png` - Differential abundance
- `Figure6.png` - Healthy vs. dysbiotic comparison
- `Figure7.png` - Fresh vs. frozen comparison
- `Figure8.png` - Anatomical site enrichment
- `Figure9.png` - Multivariate analysis

### Supplementary Figures
- `SupplementaryFigure1.png` through `SupplementaryFigure16.png`
- Various `.csv` files with statistical results and summary tables

## Data Files

### Input Data Description

| File | Description |
|------|-------------|
| `otu_table.csv` | Operational Taxonomic Unit (OTU) abundance matrix |
| `sample_metadata.csv` | Sample information: location, antibody, condition, storage |
| `taxonomy.csv` | Taxonomic classification (Kingdom to Species) |
| `facs_binding_events.csv` | Flow cytometry binding event counts |
| `anatomical_site_enrichment.csv` | Site-specific enrichment data |
| `lactobacillus_fresh_frozen.csv` | Lactobacillus abundance comparisons |
| `seq_clean.csv` | Quality-filtered sequence data |

**Data Availability**: Raw sequencing data will be deposited to SRA/ENA upon manuscript acceptance.

## Analysis Methods

### Statistical Tests

- **Alpha diversity**: Kruskal-Wallis H-test with Dunn's post-hoc
- **Beta diversity**: PERMANOVA (permutational ANOVA)
- **Differential abundance**: Mann-Whitney U test with Bonferroni correction
- **Correlation**: Spearman rank correlation

### Color Scheme

Consistent color coding is used across all figures:

- **Anatomical Location**: Vagina (blue), Cervix (green), Endometrium (red)
- **Antibody Type**: Pre-sorted (gray), IgA (orange), IgM (purple), IgG (teal)
- **Condition**: Healthy (blue), Dysbiotic (red)
- **Storage**: Fresh (orange), Frozen (blue)

## Troubleshooting

### Common Issues

**1. Module Not Found Error**
```
ModuleNotFoundError: No module named 'numpy'
```
**Solution**: Ensure virtual environment is activated and dependencies are installed:
```bash
pip install -r requirements.txt
```

**2. File Not Found Error**
```
FileNotFoundError: [Errno 2] No such file or directory: '../data/otu_table.csv'
```
**Solution**: Run scripts from the repository root directory, not from the `scripts/` subdirectory.

**3. Permission Denied (Windows)**
```
Activate.ps1 cannot be loaded because running scripts is disabled
```
**Solution**: Set PowerShell execution policy:
```powershell
Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser
```

**4. Memory Error**
**Solution**: Close other applications or use a system with more RAM.


## Contact

For questions, issues, or collaboration inquiries:

- **Prakash Lingasamy**: prakash.lingasamy@ut.ee
- **Andres Salumets**: andres.salumets@ki.se

## License

This project is intended for academic and research purposes. For commercial use or collaborations, please contact the authors.

## Acknowledgments

We thank all contributors to the mFLOW-Seq methodology development and the research participants who made this study possible.