"""
00_shared_utils.py – Shared helper functions for mFLOW-Seq bioinformatics figures
==================================================================================
Import this module from each figure script.

Functions:
  load_data()          → (otu_rel, meta, tax)  — normalised OTU + metadata + taxonomy
  bray_curtis(X)       → distance matrix (numpy)
  pcoa(D)              → coordinates, variance explained
  shannon(v)           → float
  simpson(v)           → float
  richness(v)          → int
  hellinger(X)         → transformed count matrix
  permanova(D, groups) → dict with F, p, R2
  adonis_pairwise(...) → DataFrame with pairwise PERMANOVA results
  kw_pairwise(data, groups, alpha) → DataFrame with Kruskal-Wallis + Dunn post-hoc
  pct_label(p)         → significance stars string

Colour conventions (consistent across all figures):
  LOC_C   — Location  (Vagina / Cervix / Endometrium)
  AB_C    — Antibody  (PRE / A / M / G)
  COND_C  — Condition (HEALTHY / DYSBIOTIC)
  STOR_C  — Storage   (FRESH / FROZEN)
"""

import os
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist, squareform
from scipy.stats import kruskal, mannwhitneyu

# ── Colour palettes ───────────────────────────────────────────────────────────
LOC_C  = {'VAGINA': '#2166AC', 'CERVIX': '#4DAC26', 'ENDOMETRIUM': '#D7191C'}
AB_C   = {'PRE': '#636363', 'A': '#E66101', 'M': '#5E3C99', 'G': '#1B9E77'}
COND_C = {'HEALTHY': '#2196F3', 'DYSBIOTIC': '#F44336'}
STOR_C = {'FRESH': '#FF8C00', 'FROZEN': '#1565C0'}

LOC_ORDER  = ['VAGINA', 'CERVIX', 'ENDOMETRIUM']
AB_ORDER   = ['PRE', 'A', 'M', 'G']
AB_LABELS  = {'PRE': 'Pre-sorted', 'A': 'IgA', 'M': 'IgM', 'G': 'IgG'}

# ── Paths ─────────────────────────────────────────────────────────────────────
_HERE      = os.path.dirname(os.path.abspath(__file__))
DATA_DIR   = os.path.join(_HERE, '..', 'data')
FIG_DIR    = os.path.join(_HERE, '..', 'new_figures')
os.makedirs(FIG_DIR, exist_ok=True)


# ── Data loading ──────────────────────────────────────────────────────────────
def load_data():
    """
    Load and normalise data.

    Returns
    -------
    otu_rel : pd.DataFrame  — relative abundance (0-1), shape (taxa × samples)
    meta    : pd.DataFrame  — sample metadata, indexed by SampleID
    tax     : pd.DataFrame  — taxonomy table, indexed by Taxon
    """
    # otu_table.csv: Burkholderia cepacia and Streptococcus pneumoniae
    # removed (confirmed reagent/kit contaminants). Taxa <0.5% are intentionally
    # retained — they contribute <0.4% of total reads and removal did not alter
    # community-level results; retaining them is the more transparent approach.
    otu  = pd.read_csv(os.path.join(DATA_DIR, 'otu_table.csv'), index_col=0)
    meta = pd.read_csv(os.path.join(DATA_DIR, 'sample_metadata.csv'),
                       index_col='SampleID')
    tax  = pd.read_csv(os.path.join(DATA_DIR, 'taxonomy.csv'),
                       index_col='Taxon')

    # Normalise to relative abundance
    col_sums = otu.sum(axis=0)
    col_sums[col_sums == 0] = 1
    otu_rel = otu.div(col_sums, axis=1)

    # Keep only samples present in both OTU and metadata
    common = otu_rel.columns.intersection(meta.index)
    otu_rel = otu_rel[common]
    meta    = meta.loc[common]

    return otu_rel, meta, tax


# ── Diversity metrics ─────────────────────────────────────────────────────────
def shannon(v):
    """Shannon entropy H' (base e)."""
    v = np.array(v, dtype=float)
    v = v[v > 0]
    if len(v) == 0:
        return 0.0
    p = v / v.sum()
    return float(-np.sum(p * np.log(p)))


def simpson(v):
    """Simpson diversity (1 – D)."""
    v = np.array(v, dtype=float)
    v = v[v > 0]
    if len(v) == 0:
        return 0.0
    p = v / v.sum()
    return float(1 - np.sum(p ** 2))


def richness(v):
    """Observed species richness."""
    return int(np.sum(np.array(v) > 0))


def evenness(v):
    """Pielou's evenness (J' = H' / ln(S))."""
    v = np.array(v, dtype=float)
    s = richness(v)
    h = shannon(v)
    if s <= 1:
        return 0.0
    return float(h / np.log(s))


# ── Transformations ───────────────────────────────────────────────────────────
def hellinger(X):
    """
    Hellinger transformation: sqrt(relative abundance).
    X : array-like, shape (n_samples, n_taxa)
    """
    X = np.array(X, dtype=float)
    row_sums = X.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1
    return np.sqrt(X / row_sums)


def clr(X, pseudocount=0.5):
    """
    Centred log-ratio transformation.
    X : array-like, shape (n_samples, n_taxa)
    """
    X = np.array(X, dtype=float) + pseudocount
    log_X = np.log(X)
    return log_X - log_X.mean(axis=1, keepdims=True)


# ── Beta diversity ─────────────────────────────────────────────────────────────
def bray_curtis(X):
    """
    Bray-Curtis dissimilarity matrix.
    X : array-like, shape (n_samples, n_taxa) — relative abundances.
    Returns square distance matrix as numpy array.
    """
    X  = np.array(X, dtype=float)
    n  = X.shape[0]
    D  = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            num = np.sum(np.abs(X[i] - X[j]))
            den = np.sum(X[i]) + np.sum(X[j])
            d   = num / den if den > 0 else 0.0
            D[i, j] = D[j, i] = d
    return D


def pcoa(D, n_axes=4):
    """
    Principal Coordinates Analysis (PCoA).

    Parameters
    ----------
    D      : square distance matrix (numpy array)
    n_axes : number of axes to return

    Returns
    -------
    coords       : numpy array, shape (n, n_axes)
    explained_var: list of variance explained per axis (0-1)
    """
    n   = D.shape[0]
    D2  = D ** 2
    J   = np.eye(n) - np.ones((n, n)) / n
    B   = -0.5 * J @ D2 @ J
    B   = (B + B.T) / 2          # enforce symmetry

    eigvals, eigvecs = np.linalg.eigh(B)
    order   = np.argsort(eigvals)[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]

    # Keep only positive eigenvalues
    pos     = eigvals > 1e-10
    coords  = np.zeros((n, n_axes))
    var_exp = []
    total   = eigvals[pos].sum() if pos.any() else 1.0

    for k in range(n_axes):
        if k < pos.sum():
            coords[:, k] = eigvecs[:, k] * np.sqrt(eigvals[k])
            var_exp.append(eigvals[k] / total)
        else:
            var_exp.append(0.0)

    return coords, var_exp


def nmds(D, n_axes=2, n_init=10, max_iter=300, seed=42):
    """
    Non-metric MDS (wrapper around sklearn).
    D : square distance matrix.
    """
    from sklearn.manifold import MDS
    mds = MDS(n_components=n_axes, metric=False, dissimilarity='precomputed',
              n_init=n_init, max_iter=max_iter, random_state=seed, n_jobs=1)
    return mds.fit_transform(D)


# ── Statistical tests ─────────────────────────────────────────────────────────
def permanova(D, groups, n_perm=999, seed=42):
    """
    PERMANOVA (adonis) for a single factor.

    Parameters
    ----------
    D      : square distance matrix (numpy)
    groups : array-like of group labels, same order as D rows/cols
    n_perm : number of permutations

    Returns
    -------
    dict with keys: F, R2, p, df_model, df_resid
    """
    np.random.seed(seed)
    groups = np.array(groups)
    n      = len(groups)
    unique = np.unique(groups)
    k      = len(unique)

    def _ss_total(D):
        return (D ** 2).sum() / (2 * n)

    def _ss_within(D, grps):
        ss = 0.0
        for g in np.unique(grps):
            idx = np.where(grps == g)[0]
            sub = D[np.ix_(idx, idx)]
            ni  = len(idx)
            if ni > 1:
                ss += (sub ** 2).sum() / (2 * ni)
        return ss

    ss_tot = _ss_total(D)
    ss_w   = _ss_within(D, groups)
    ss_b   = ss_tot - ss_w
    df_b   = k - 1
    df_w   = n - k
    F_obs  = (ss_b / df_b) / (ss_w / df_w) if df_w > 0 else np.nan
    R2_obs = ss_b / ss_tot if ss_tot > 0 else np.nan

    # Permutation
    count = 0
    for _ in range(n_perm):
        perm_grps = np.random.permutation(groups)
        ss_w_p    = _ss_within(D, perm_grps)
        ss_b_p    = ss_tot - ss_w_p
        F_p       = (ss_b_p / df_b) / (ss_w_p / df_w) if df_w > 0 else np.nan
        if not np.isnan(F_p) and F_p >= F_obs:
            count += 1

    p_val = (count + 1) / (n_perm + 1)
    return {'F': F_obs, 'R2': R2_obs, 'p': p_val,
            'df_model': df_b, 'df_resid': df_w}


def pct_label(p):
    """Return significance asterisks for p-value."""
    if p < 0.001: return '***'
    if p < 0.01:  return '**'
    if p < 0.05:  return '*'
    return 'ns'
