# mFLOW-Seq revision v25.2.1

**Immunoglobulin class–resolved engagement of membrane-intact microbiota across female reproductive tract compartments using mFLOW-Seq**

Analysis scripts, original coded inputs, saved resampling records, numerical results and selected figures accompanying the Microbiome revision. The retained sequencing dataset comprises 149 libraries from 15 women and 40 recorded preparations. Participant, preparation and library are different units.

## Start here

- [Reproduce and verify](REPRODUCE.md): portable integrity checks and the existing analysis/figure entry points.
- [Current supplementary workbook](revisions/revision_v25_2_1/selected_tables/Supplementary_Tables_CLEAN.xlsx): the exact 39-sheet CLEAN workbook from v25.2.1.
- [Current figures and captions](revisions/revision_v25_2_1/selected_figures): six main and nineteen numbered supplementary figures, including multi-part control displays; original microscopy bytes retained.
- [Publication source index](revisions/revision_v25_2_1/Publication_table_source_index.csv): 82 source blocks, extending the original 78 with four supplied sequencing-QC/benchmark tables. Paths and SHA256 values identify the preserved source CSVs.
- [Sequencing inputs and commands](revisions/revision_v25): supplied 149-library accounting workbook/CSV, three computational benchmark types, main.nf, run.sh and arithmetic checks.
- [Existing downstream analysis sources](revisions/revision_v24_8/source_materials): participant/preparation-aware models, site-specific scores, Q15 joint model and own-PRE comparisons, FACS source counts and figure-reproduction scripts.

The v25.2.1 manuscript polish did not change scientific data, analyses or figure pixels. Older source directories retain their historical filenames and README text; use this root README and REPRODUCE.md for the current entry point. The earlier root README and its matching checksum manifest are preserved under historical_documentation. The original published repository's older run_all.py/data/Output content remains historical and is not the entry point for the revised figures. Apply this candidate as an addition/update to the existing repository; do not delete its older code or data.

## Sequencing and evidence scope

The confirmed workflow uses 27F/1492R near-full-length V1–V9 amplification, Tn5-fragmented libraries and NextSeq 1000 single-end reads delivered after service-side adapter/quality processing. The deposited reads are 80 nt; the Bracken -r 75 setting selects a reference distribution and is not the delivered read length. The analyst identifies Kraken2 Standard k2_standard_20221209 and AbundanceMi/internal-reference version 2.3. Supplied workflow commands and summaries improve transparency but do not reconstruct the complete original proprietary FASTQ-to-final-matrix processing chain.

Public SRA access under PRJNA1454428/SRP692419 was checked on 14 November 2026: 149 runs and representative sequence content were accessible. Deposited FASTQs are provider-preprocessed, not unprocessed instrument output. A complete transfer/checksum verification was not performed. Older statements that the reads remain embargoed are superseded. Minor public metadata maintenance is separate from this local release preparation.

Physical D6310 observations are three independent extraction replicates sharing downstream batches. The single/equal/pair sequence simulations are separate computational benchmarks. Species calls remain provisional; mixture estimates show substantial distortion. Neither a benchmark table nor successful downstream replay establishes uniformly accurate clinical species identification.

DAPI-negative sorting is operational selection intended to enrich for putatively membrane-intact bacteria. HUT093 BSA/PBS1/PBS2 contain DAPI, no antibodies and no biological specimen; their fluorescence-region occupancies are background observations, not IgM-bound bacteria or clinical false-positive rates. Retrospective reports and historical statistics remain distinct. These data do not establish quantitative clinical DAPI agreement, matched antibody specificity, viability, contamination-free sampling, biological ascension or cryopreservation equivalence.

## Version and access

Prepared version: **v25.2.1**. Existing repository: https://github.com/Prakashbio/mflow-seq-immune-microbiome-reproductive-tract . Existing Zenodo record supplied with the study: https://doi.org/10.5281/zenodo.19437248 . These links identify existing public materials; this candidate has not been uploaded or assigned a new version DOI. After authorized publication, record the actual release tag/DOI without inventing one in advance.

## License and provenance

The existing GPLv3 code LICENSE is retained unchanged. Original code, source files, data licenses and attribution are preserved; no new blanket data license is assigned by this preparation. Original source copies are not silently rewritten to match later annotations. The preservation ledger records every member of the previous code archive and the two historical-document relocations. New journal-facing tables and figure selection are clearly separated from those older source files.
