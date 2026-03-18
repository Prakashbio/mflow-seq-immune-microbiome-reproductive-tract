# mFLOW-Seq Pipeline for Immune-Coated Microbiome Analysis in the Female Reproductive Tract

[![Python](https://img.shields.io/badge/Python-3.9+-blue.svg)](https://www.python.org/)
[![License: GPL v3](https://img.shields.io/badge/License-Research-orange.svg)](LICENSE)

Bioinformatics analysis pipeline and figure generation scripts for the mFLOW-Seq microbiome study investigating spatial immunoglobulin engagement of live microbiota across female reproductive tract compartments.

This repository contains all scripts required to reproduce the main and supplementary analyses presented in the manuscript. All figures can be regenerated directly from the provided code and processed datasets.

## Overview

This repository contains the complete computational analysis workflow used to generate all main and supplementary figures for the manuscript:

> **Spatial Immunoglobulin Engagement of Live Microbiota Across Female Reproductive Tract Compartments Revealed by mFLOW-Seq**
>
> Prakash Lingasamy, Merli Saare, Kasun Madhuranga Godakuma Godagedara, Sergio Vela Moreno, Karita Särekannu, Dmitri Lubenets, Alberto Sola-Leyva, Naisarg Patel, Vijayachitra Modhukur, Reet Mändar, Andres Salumets
>
> *Manuscript under communication*

### About mFLOW-Seq

mFLOW-Seq (Microbiome Flow-Sorting and Sequencing) is a novel methodology combining fluorescence-activated cell sorting (FACS) and 16S rRNA sequencing to profile immunoglobulin-bound bacterial populations in complex microbiome samples. This approach enables the characterization of host-microbiota interactions t the level of FACS-sorted bacterial subpopulations.

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

The pipeline generates all main and supplementary figures from the manuscript. Below is a complete description of every script and the analyses it performs.

### Main Figures

| Script | Figure | Description | Key analyses |
|--------|--------|-------------|--------------|
| `fig2_alpha_diversity.py` | Figure 2 | Alpha diversity across Ig fractions and anatomical sites, stratified by clinical condition | Shannon H′, Simpson (1−D), Pielou's evenness, Kruskal–Wallis, Mann–Whitney U |
| `fig3_beta_diversity.py` | Figure 3 | Beta diversity ordination: PCoA and NMDS of Bray–Curtis dissimilarities by location, Ig fraction, condition, and storage | Bray–Curtis, Hellinger transform, PCoA, NMDS, PERMANOVA (adonis2) |
| `fig4_facs_binding.py` | Figure 4 | FACS-sorted bacterial event counts and 16S rRNA sequencing read distributions per Ig fraction and anatomical site | Violin plots, Wilcoxon signed-rank test (IgA vs IgM) |
| `fig5_differential_abundance.py` | Figure 5 | Differential Ig-coating of the microbiota: heatmap, log₂(IgA/IgM) diverging bars, and IgA vs IgM bubble chart | IgA coating score, log₂ fold-change visualisation, taxon-level enrichment heatmap |
| `fig6_healthy_vs_dysbiotic.py` | Figure 6 | Healthy vs. microbiologically dysbiotic community comparison: key taxon abundances, Shannon diversity, PCoA | Mann–Whitney U, Bray–Curtis PCoA, PERMANOVA per site |
| `fig7_fresh_vs_frozen.py` | Figure 7 | Fresh vs. cryopreserved sample comparison: *Lactobacillus* spp. abundance and taxon composition heatmaps | Mann–Whitney U, taxon abundance heatmaps |
| `fig8_anatomical_site.py` | Figure 8 | IgA vs. IgM taxon enrichment lollipop plots across vaginal, cervical, and endometrial compartments | IgA − IgM raw abundance difference, contaminant flagging |
| `fig9_multivariate_analysis.py` | Figure 9 | Advanced multivariate analysis: CLR-PCA biplot, PERMANOVA effect-size barplot, and spatial Ig gradient for key taxa | CLR transformation, PCA biplot, PERMANOVA R² effect sizes, bar chart summaries |

### Supplementary Figures

| Script | Figure | Description | Key analyses |
|--------|--------|-------------|--------------|
| `sfig1_phylum_composition.py` | S1 | Phylum-level bacterial composition across Ig fractions, anatomical sites, and clinical condition | Stacked bar plots, relative abundance |
| `sfig2_alpha_diversity_storage.py` | S2 | Alpha diversity stratified by storage method (fresh vs. frozen) across all Ig fractions and sites | Shannon H′, Simpson, Pielou's evenness, Kruskal–Wallis, Mann–Whitney U |
| `sfig3_per_location_pcoa.py` | S3 | Per-anatomical-site Bray–Curtis PCoA coloured by Ig fraction | Bray–Curtis, PCoA, PERMANOVA per site |
| `sfig4_vs_presorted.py` | S4 | Compositional comparison of Ig-sorted fractions vs. pre-sorted baseline (log₂ enrichment) | log₂(Ig fraction/Pre-sorted) per taxon, enrichment plots |
| `sfig5_condition_volcano.py` | S5 | Volcano plots of healthy vs. dysbiotic differential taxa per anatomical site | Mann–Whitney U, −log₁₀(p) vs log₂ fold-change |
| `sfig6_cooccurrence_matrix.py` | S6 | Taxon co-occurrence correlation matrices per anatomical site | Spearman rank correlation, hierarchical clustering heatmaps |
| `sfig7_fresh_frozen_correlation.py` | S7 | Fresh vs. frozen taxon abundance scatter plots with Spearman r per site | Spearman rank correlation, mean relative abundance |
| `sfig8_taxon_accumulation.py` | S8 | Taxon accumulation (rarefaction) curves by anatomical site, confirming sampling completeness | Random rarefaction subsampling, 95% confidence intervals |
| `sfig9_supplementary_igg.py` | S9 | IgG enrichment lollipop plots, proportional Ig-coating space (IgA/IgM/IgG), and pairwise Bray–Curtis distances between fractions | log₂(IgG/Pre-sorted), proportional coating ratios, Bray–Curtis within-fraction distances |
| `sfig10_iga_coating_score.py` | S10 | Per-taxon IgA coating scores [IgA/(IgA+IgM)] across anatomical sites, *Lactobacillus* vs. other genera, and healthy vs. dysbiotic comparisons | IgA coating score (formula-based), dot-plot heatmap, Mann–Whitney U |
| `sfig11_patient_paired_readdepth.py` | S11 | Per-patient alpha diversity trajectories across sites, read depth QC, within- vs between-patient community similarity, read depth vs Shannon correlation | Spearman r (depth vs diversity), Bray–Curtis pairwise distances, Mann–Whitney U |
| `sfig12_lactobacillus_gradient.py` | S12 | *Lactobacillus* species composition across sites and fractions, species-level enrichment vs pre-sorted, and dysbiotic taxon IgA/IgM abundance by condition | log₂(fraction/pre-sorted), bar plots by species and condition |
| `sfig13_composition.py` | S13 | Species-level microbial community composition (top 20 taxa) across all Ig fractions, sites, and clinical conditions | Stacked bar plots, mean relative abundance |
| `sfig14_prevalence_network.py` | S14 | Filtering threshold comparison, taxon prevalence heatmap, and microbial co-occurrence network (pre-sorted) | Taxon prevalence scoring, Spearman correlation network (\|r\| ≥ 0.4, p < 0.05) |
| `sfig15_fresh_frozen.py` | S15 | Community structure and alpha diversity stratified by storage condition per anatomical site | Bray–Curtis PCoA, PERMANOVA, Shannon H′, Mann–Whitney U |
| `sfig16_contamination_analysis.py` | S16 | Contamination analysis: decontam-equivalent prevalence scoring and Bayesian source tracking (SourceTracker framework) | Prevalence-based decontam scoring, Mann–Whitney U filter validation |

### Shared Utilities (`shared_utils.py`)

Core functions used by all figure scripts:

| Function | Description |
|----------|-------------|
| `load_data()` | Loads decontaminated OTU table, metadata, and taxonomy |
| `shannon(v)` | Shannon entropy H′ = −Σ pᵢ ln pᵢ |
| `simpson(v)` | Simpson's diversity 1−D |
| `richness(v)` | Observed taxon richness |
| `evenness(v)` | Pielou's evenness J′ = H′/ln(S) |
| `hellinger(X)` | Hellinger transformation: √(relative abundance) |
| `clr(X)` | Centred log-ratio transformation |
| `bray_curtis(X)` | Bray–Curtis dissimilarity matrix |
| `pcoa(D)` | Principal Coordinates Analysis |
| `nmds(D)` | Non-metric Multidimensional Scaling (sklearn wrapper) |
| `permanova(D, groups)` | Pseudo-F PERMANOVA with permutation p-value |
| `pct_label(p)` | Significance stars (\*, \*\*, \*\*\*, ns) |

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
git clone https://github.com/Prakashbio/mflow-seq-immune-microbiome-reproductive-tract
cd mflow-seq-immune-microbiome-reproductive-tract
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

**Note**: `run_all.py` executes each script with its working directory set to `scripts/`. Individual scripts can also be run directly from the `scripts/` directory.

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
| `otu_table.csv` | **Primary analysis file** — decontaminated OTU abundance matrix (39 taxa, 149 samples; *Burkholderia cepacia* and *Streptococcus pneumoniae* removed) |
| `sample_metadata.csv` | Sample information: location, antibody, condition, storage |
| `taxonomy.csv` | Taxonomic classification (Kingdom to Species) |
| `facs_binding_events.csv` | Flow cytometry binding event counts |
| `anatomical_site_enrichment.csv` | Site-specific enrichment data |
| `lactobacillus_fresh_frozen.csv` | Lactobacillus abundance comparisons |
| `seq_clean.csv` | Quality-filtered sequence data |
| `otu_clean_decontaminated.csv` | Same as `otu_table.csv` — provided for cross-reference |
| `otu_clean.csv` | Pre-decontamination OTU table (41 taxa; includes *Burkholderia cepacia* and *Streptococcus pneumoniae* before removal) |
| `taxa_abundance_fresh_frozen.csv` | Mean taxon abundance per storage condition and site |
| `filtering_statistics.xlsx` | Taxon counts at different abundance thresholds per group |

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

This project is licensed under the GNU General Public License v3.0 (GPL-3.0). See the [LICENSE](LICENSE) file for details.

If you use this code in your research, please cite the associated manuscript (see [Citation](#citation) below).

## Citation

If you use this code or data, please cite:

> Lingasamy P, Saare M, Godagedara KMG, Vela Moreno S, Särekannu K, Lubenets D, Sola-Leyva A, Modhukur V, Mändar R, Salumets A.
> **Spatial Immunoglobulin Engagement of Live Microbiota Across Female Reproductive Tract Compartments Revealed by mFLOW-Seq.**
> *Manuscript under review.*

## Acknowledgments

We thank all contributors to the mFLOW-Seq methodology development and the research participants who made this study possible.