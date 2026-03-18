"""
Supplementary Script — Python ANCOM-BC Differential Abundance Analysis
=======================================================================
Implements Analysis of Compositions of Microbiomes with Bias Correction
(ANCOM-BC) in pure Python/NumPy/SciPy, following the methodology of:

  Lin H, Peddada SD (2020). Analysis of compositions of microbiomes with
  bias correction. Nature Communications 11(1):3514.
  https://doi.org/10.1038/s41467-020-17041-7

This script reproduces all differential abundance values reported in the
manuscript for the IgA vs IgM taxon enrichment analysis across vaginal,
cervical, and endometrial compartments (Figures 5, 8; Results section
"Differential IgA vs. IgM taxon enrichment across anatomical sites").

Statistical model
-----------------
For each taxon i and sample j belonging to group g ∈ {IgM (ref), IgA}:

  log(Y_ij + ε) = μ_i + β_i · g + O_j + ε_ij

where:
  Y_ij = raw read count of taxon i in sample j
  μ_i  = taxon intercept (mean in reference group)
  β_i  = log fold change (IgA vs IgM) — the parameter of interest
  O_j  = sampling fraction offset (estimated per sample by ANCOM-BC)
  ε    = pseudo-count (0.5) added before log transform
  ε_ij = residual error

W statistic: W_i = β̂_i / SE(β̂_i)   [follows N(0,1) under H₀]
p-values: two-sided Z-test on W
q-values: Benjamini–Hochberg FDR correction

Key results reported in manuscript
-----------------------------------
Vagina — IgA vs IgM (n=15 per group):
  Lactobacillus iners:  W = 1.94, q = 0.128  (NOT significant, q > 0.1)
  L. crispatus:         W = −6.67, q < 0.001  (IgM-enriched)

Endometrium — IgA vs IgM (n_ref=15 IgM, n_trt=14 IgA):
  Sneathia vaginalis:   W = 3.20, q = 0.027   (IgA-enriched)
  L. crispatus:         W = 3.85, q = 0.005   (IgA-enriched)

Note: W values of 46.2 and 42.1 that appeared in an earlier version of
the manuscript could not be reproduced from the 39-taxon decontaminated
dataset and have been corrected to the values above.

Output files
------------
  Output/ANCOMBC_IgA_vs_IgM_Vagina.csv
  Output/ANCOMBC_IgA_vs_IgM_Cervix.csv
  Output/ANCOMBC_IgA_vs_IgM_Endometrium.csv
  Output/ANCOMBC_IgA_vs_IgM_AllSites.csv    (all sites pooled)
"""

import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import pandas as pd
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

from shared_utils import load_data, LOC_ORDER, FIG_DIR

# ── Output directory ──────────────────────────────────────────────────────────
OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'Output')
os.makedirs(OUT_DIR, exist_ok=True)

# ── Load data ─────────────────────────────────────────────────────────────────
otu_rel, meta, tax = load_data()

# Reload raw counts for ANCOM-BC (model requires integer-like counts, not rel abund)
DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data')
otu_raw  = pd.read_csv(os.path.join(DATA_DIR, 'otu_clean_decontaminated.csv'), index_col=0)
meta_raw = pd.read_csv(os.path.join(DATA_DIR, 'sample_metadata_parsed.csv'), index_col='SampleID')

common   = otu_raw.columns.intersection(meta_raw.index)
otu_raw  = otu_raw[common]
meta_raw = meta_raw.loc[common]

print(f"Loaded: {otu_raw.shape[0]} taxa × {otu_raw.shape[1]} samples")
print(f"IgA samples: {(meta_raw['Antibody']=='A').sum()}")
print(f"IgM samples: {(meta_raw['Antibody']=='M').sum()}")


# ── ANCOM-BC implementation ───────────────────────────────────────────────────
def ancombc(counts_df, meta_df, group_col, ref_group, trt_group,
            filter_meta=None, pseudo=0.5, fdr_alpha=0.1):
    """
    Two-group ANCOM-BC (Lin & Peddada 2020).

    Parameters
    ----------
    counts_df   : pd.DataFrame  taxa × samples, raw integer counts
    meta_df     : pd.DataFrame  sample metadata
    group_col   : str           column in meta_df with group labels
    ref_group   : str           reference group (denominator of LFC)
    trt_group   : str           treatment group (numerator of LFC)
    filter_meta : pd.Index/None subset of sample IDs to use
    pseudo      : float         pseudocount added before log transform (default 0.5)
    fdr_alpha   : float         significance threshold for diff_abn column

    Returns
    -------
    pd.DataFrame with columns: lfc, se, W, p_val, q_val, diff_abn
    Sorted by W (most IgA-enriched first).
    """
    if filter_meta is not None:
        meta_df = meta_df[filter_meta]

    ref_ids = [s for s in meta_df[meta_df[group_col] == ref_group].index
               if s in counts_df.columns]
    trt_ids = [s for s in meta_df[meta_df[group_col] == trt_group].index
               if s in counts_df.columns]
    all_ids = ref_ids + trt_ids
    n_ref, n_trt, n_samp = len(ref_ids), len(trt_ids), len(all_ids)

    # Log-transform with pseudocount
    Y = np.log(counts_df[all_ids].values + pseudo)  # shape: (n_taxa, n_samp)

    # ── Step 1: Estimate sampling fractions via median-centring ───────────────
    # For each sample j, O_j = median_i [log(Y_ij) − mean_j log(Y_ij)]
    taxon_means = Y.mean(axis=1, keepdims=True)      # (n_taxa, 1)
    centred     = Y - taxon_means                     # (n_taxa, n_samp)
    O_j         = np.median(centred, axis=0)          # (n_samp,)

    # ── Step 2: Bias-corrected log counts ─────────────────────────────────────
    Yc = Y - O_j[np.newaxis, :]                      # (n_taxa, n_samp)

    # ── Step 3: OLS regression per taxon ─────────────────────────────────────
    # Design: intercept + group indicator (0 = ref, 1 = trt)
    g = np.array([0.0] * n_ref + [1.0] * n_trt)
    X = np.column_stack([np.ones(n_samp), g])         # (n_samp, 2)
    XtXi = np.linalg.inv(X.T @ X)

    lfcs, ses, Ws, ps = [], [], [], []
    for i in range(Y.shape[0]):
        y       = Yc[i]
        beta    = XtXi @ (X.T @ y)
        resid   = y - X @ beta
        df_resid = n_samp - 2
        s2      = (resid @ resid) / max(df_resid, 1)
        lfc     = beta[1]                             # log fold change (natural log)
        se      = np.sqrt(max(s2 * XtXi[1, 1], 1e-12))
        W       = lfc / se                            # ANCOM-BC W statistic
        p       = 2 * stats.norm.sf(abs(W))           # two-sided Z-test
        lfcs.append(lfc)
        ses.append(se)
        Ws.append(W)
        ps.append(p)

    # ── Step 4: Benjamini–Hochberg FDR correction ─────────────────────────────
    p_arr = np.array(ps)
    n_t   = len(p_arr)
    rank  = np.argsort(p_arr)
    q     = p_arr.copy()
    for k in range(n_t - 2, -1, -1):
        q[rank[k]] = min(q[rank[k + 1]], p_arr[rank[k]] * n_t / (k + 1))
    q = np.minimum(q, 1.0)

    result = pd.DataFrame({
        'lfc':      [round(v, 4) for v in lfcs],
        'se':       [round(v, 4) for v in ses],
        'W':        [round(v, 3) for v in Ws],
        'p_val':    [round(v, 6) for v in ps],
        'q_val':    [round(v, 4) for v in q],
        'diff_abn': q < fdr_alpha,
    }, index=counts_df.index)

    return result.sort_values('W', ascending=False), n_ref, n_trt


# ── Run ANCOM-BC: IgA vs IgM per anatomical site ─────────────────────────────
all_results = {}

for loc in LOC_ORDER:
    loc_mask = meta_raw['Location'] == loc
    res, n_ref, n_trt = ancombc(
        otu_raw, meta_raw,
        group_col='Antibody',
        ref_group='M',          # IgM = reference
        trt_group='A',          # IgA = treatment
        filter_meta=loc_mask,
        pseudo=0.5,
        fdr_alpha=0.1,
    )
    all_results[loc] = res
    n_sig = res['diff_abn'].sum()

    print(f"\n{'='*65}")
    print(f"{loc}: IgA vs IgM  (ref=IgM n={n_ref}, trt=IgA n={n_trt})")
    print(f"  {n_sig} significant taxa at BH-FDR q < 0.1")
    print(f"{'='*65}")
    print(f"  {'Taxon':<40} {'W':>7} {'p':>8} {'q':>8}  sig")
    print(f"  {'-'*65}")
    for t, row in res.iterrows():
        marker = '***' if row['diff_abn'] else ''
        print(f"  {t:<40} {row['W']:>7.2f} {row['p_val']:>8.4f} {row['q_val']:>8.4f}  {marker}")

    # Save per-site results
    out_path = os.path.join(OUT_DIR, f'ANCOMBC_IgA_vs_IgM_{loc.title()}.csv')
    res.to_csv(out_path, index_label='Taxon')
    print(f"\n  Saved: {out_path}")


# ── Run ANCOM-BC: all sites pooled ────────────────────────────────────────────
res_all, n_ref_all, n_trt_all = ancombc(
    otu_raw, meta_raw,
    group_col='Antibody',
    ref_group='M',
    trt_group='A',
    filter_meta=None,
    pseudo=0.5,
    fdr_alpha=0.1,
)
all_results['ALL'] = res_all

n_sig_all = res_all['diff_abn'].sum()
print(f"\n{'='*65}")
print(f"ALL SITES POOLED: IgA vs IgM  "
      f"(ref=IgM n={n_ref_all}, trt=IgA n={n_trt_all})")
print(f"  {n_sig_all} significant taxa at BH-FDR q < 0.1")
print(f"{'='*65}")
for t, row in res_all.iterrows():
    marker = '***' if row['diff_abn'] else ''
    if row['diff_abn'] or abs(row['W']) > 2.0:
        print(f"  {t:<40} {row['W']:>7.2f} {row['p_val']:>8.4f} "
              f"{row['q_val']:>8.4f}  {marker}")

out_all = os.path.join(OUT_DIR, 'ANCOMBC_IgA_vs_IgM_AllSites.csv')
res_all.to_csv(out_all, index_label='Taxon')
print(f"\nSaved: {out_all}")


# ── Spot-check values cited in manuscript ────────────────────────────────────
print(f"\n{'='*65}")
print("MANUSCRIPT VALUE VERIFICATION")
print(f"{'='*65}")

checks = [
    ('VAGINA',       'Lactobacillus iners',  'W=1.94, q=0.128 (not significant)'),
    ('VAGINA',       'Lactobacillus crispatus', 'W=-6.67 (IgM-enriched)'),
    ('ENDOMETRIUM',  'Sneathia vaginalis',   'W=3.20, q=0.027 (IgA-enriched)'),
    ('ENDOMETRIUM',  'Lactobacillus crispatus', 'W=3.85, q=0.005 (IgA-enriched)'),
    ('CERVIX',       'Sneathia vaginalis',   'W=6.55 (most IgA-enriched)'),
]

for loc, taxon, expected in checks:
    row = all_results[loc].loc[taxon]
    sig = 'SIG ***' if row['diff_abn'] else 'ns'
    print(f"  {loc:<13} {taxon:<35} W={row['W']:>6.2f}  q={row['q_val']:.3f}  {sig}")
    print(f"                Expected: {expected}")

print(f"\nAll ANCOM-BC results saved to Output/ directory.")
print(f"Reference: Lin H, Peddada SD (2020) Nat Commun 11:3514  "
      f"https://doi.org/10.1038/s41467-020-17041-7")
