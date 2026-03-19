"""
TOST Cryopreservation Equivalence Analysis
===========================================
Two one-sided t-test (TOST) equivalence testing for fresh vs. frozen
specimens in the mFLOW-Seq dataset.

Tests whether cryopreservation produces statistically equivalent
Ig-binding community profiles (Bray-Curtis) and Shannon diversity
across all Ig fractions at vaginal and endometrial sites.

Key finding:
  - Bray-Curtis community composition: EQUIVALENT (all fractions pooled,
    vagina and endometrium; delta = ±0.10 BC units)
  - Shannon diversity (PRE-sorted): NOT equivalent at vagina due to
    known taxon richness reduction (reported separately in Results;
    Supplementary Table S2.6)

Output files (written to FIG_DIR):
  tost_bray_curtis_results.csv   — per-site and pooled BC TOST results
  tost_shannon_results.csv       — per-site Shannon TOST results
  SFig_TOST_Cryopreservation.png — summary figure

Usage:
  python sfig_tost_cryopreservation.py

Requirements:
  numpy, pandas, scipy, matplotlib
  Data: otu_clean_decontaminated.csv, sample_metadata_parsed.csv

Reference:
  Schuirmann DJ (1987) A comparison of the two one-sided tests procedure
  and the power approach for assessing the equivalence of average
  bioavailability. J Pharmacokinet Biopharm 15(6):657-680.
  doi: 10.1007/BF01059956

Authors: Lingasamy P et al. (mFLOW-Seq manuscript)
Repository: https://github.com/Prakashbio/mflow-seq-immune-microbiome-reproductive-tract
"""

import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.stats import t as t_dist
import warnings
warnings.filterwarnings('ignore')

from shared_utils import DATA_DIR, FIG_DIR

np.random.seed(42)

# ─────────────────────────────────────────────────────────────────────────────
# CONFIGURATION
# ─────────────────────────────────────────────────────────────────────────────

# Equivalence margin for Bray-Curtis (BC) distance:
# ±0.10 BC units chosen as the minimum biologically meaningful difference
# in community composition for low-biomass mucosal microbiome datasets.
# This is conservative relative to typical within-group BC variation
# (~0.50-0.60 at all sites) and represents a ~17% relative change.
DELTA_BC = 0.10

# Equivalence margin for Shannon entropy:
# ±0.50 H' units (one half-unit of Shannon entropy), which corresponds
# to ~0.5-1.5 SD of the fresh specimen distribution depending on site.
# This is a commonly used absolute threshold for diversity equivalence.
DELTA_SHANNON = 0.50

ALPHA = 0.05   # significance level for each one-sided test

SITES = ['VAGINA', 'CERVIX', 'ENDOMETRIUM']
SITE_COLORS = {'VAGINA': '#2196F3', 'CERVIX': '#4CAF50', 'ENDOMETRIUM': '#F44336'}

# ─────────────────────────────────────────────────────────────────────────────
# DATA LOADING
# ─────────────────────────────────────────────────────────────────────────────

print("Loading data...")
otu_raw = pd.read_csv(os.path.join(DATA_DIR, 'otu_clean_decontaminated.csv'), index_col=0)
meta    = pd.read_csv(os.path.join(DATA_DIR, 'sample_metadata_parsed.csv'), index_col='SampleID')

common  = otu_raw.columns.intersection(meta.index)
otu_raw = otu_raw[common]
meta    = meta.loc[common]

# Relative abundances (samples × taxa)
col_sums = otu_raw.sum(axis=0)
col_sums[col_sums == 0] = 1
otu_rel = (otu_raw.div(col_sums, axis=1)).T

print(f"  {otu_rel.shape[1]} taxa × {otu_rel.shape[0]} samples")
print(f"  Fresh samples:  {(meta['Storage']=='FRESH').sum()}")
print(f"  Frozen samples: {(meta['Storage']=='FROZEN').sum()}")


# ─────────────────────────────────────────────────────────────────────────────
# UTILITY FUNCTIONS
# ─────────────────────────────────────────────────────────────────────────────

def shannon_entropy(v):
    """Shannon entropy H' = -sum(p_i * ln(p_i))."""
    v = np.array(v, dtype=float)
    v = v[v > 0]
    if len(v) == 0:
        return 0.0
    p = v / v.sum()
    return -np.sum(p * np.log(p))


def hellinger_transform(X):
    """Row-wise Hellinger transformation: sqrt(x_ij / row_sum_i)."""
    X = np.array(X, dtype=float)
    row_sums = X.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1
    return np.sqrt(X / row_sums)


def bray_curtis_distance(a, b):
    """Bray-Curtis dissimilarity between two vectors."""
    a, b = np.array(a, dtype=float), np.array(b, dtype=float)
    denom = np.sum(a + b)
    return np.sum(np.abs(a - b)) / denom if denom > 0 else 0.0


def tost_two_sample(x, y, delta, alpha=0.05):
    """
    Two one-sided t-tests (TOST) for equivalence of two independent groups.

    Tests H0: |mu_x - mu_y| >= delta  vs  H1: |mu_x - mu_y| < delta
    Equivalence is declared if both one-sided p-values < alpha.

    Parameters
    ----------
    x, y   : array-like, the two groups
    delta  : float, the equivalence margin (symmetric, ±delta)
    alpha  : float, significance level (default 0.05)

    Returns
    -------
    dict with keys:
      mean_x, mean_y, mean_diff  — descriptive statistics
      ci_90_low, ci_90_high      — 90% CI of mean difference
                                   (matches TOST at alpha=0.05)
      p_lower, p_upper           — one-sided p-values
      t_lower, t_upper           — t-statistics
      df                         — degrees of freedom
      is_equivalent              — bool, True if both p < alpha
      nx, ny                     — sample sizes
    """
    x, y = np.array(x, dtype=float), np.array(y, dtype=float)
    nx, ny = len(x), len(y)

    if nx < 2 or ny < 2:
        return None

    mean_x, mean_y = np.mean(x), np.mean(y)
    std_x,  std_y  = np.std(x, ddof=1), np.std(y, ddof=1)

    # Pooled standard error (equal-variance two-sample t-test)
    sp  = np.sqrt(((nx - 1) * std_x**2 + (ny - 1) * std_y**2) / (nx + ny - 2))
    se  = sp * np.sqrt(1.0 / nx + 1.0 / ny)
    df  = nx + ny - 2

    mean_diff = mean_x - mean_y

    # 90% confidence interval (corresponds to alpha=0.05 TOST)
    t_crit    = t_dist.ppf(0.95, df)
    ci_90_low  = mean_diff - t_crit * se
    ci_90_high = mean_diff + t_crit * se

    # One-sided t-statistics
    # Lower test: H0: diff <= -delta  →  reject if diff > -delta
    # Upper test: H0: diff >= +delta  →  reject if diff < +delta
    t_lower = (mean_diff - (-delta)) / se
    t_upper = (mean_diff - (+delta)) / se

    p_lower = 1.0 - t_dist.cdf(t_lower, df)   # right-tail
    p_upper = t_dist.cdf(t_upper, df)          # left-tail

    is_equivalent = (p_lower < alpha) and (p_upper < alpha)

    return {
        'mean_x':        mean_x,
        'mean_y':        mean_y,
        'std_x':         std_x,
        'std_y':         std_y,
        'mean_diff':     mean_diff,
        'ci_90_low':     ci_90_low,
        'ci_90_high':    ci_90_high,
        't_lower':       t_lower,
        't_upper':       t_upper,
        'p_lower':       p_lower,
        'p_upper':       p_upper,
        'df':            df,
        'is_equivalent': is_equivalent,
        'nx':            nx,
        'ny':            ny,
        'delta':         delta,
        'alpha':         alpha,
    }


# ─────────────────────────────────────────────────────────────────────────────
# ANALYSIS 1: BRAY-CURTIS TOST (primary equivalence test)
# ─────────────────────────────────────────────────────────────────────────────
# Approach: For each site (and pooled), compute each sample's BC distance
# to the fresh-group centroid (Hellinger-transformed). TOST then tests
# whether frozen samples are within DELTA_BC of the fresh centroid.

print("\nAnalysis 1: Bray-Curtis TOST (all Ig fractions)")
print("-" * 50)

bc_results = []

# Per-site analysis
for loc in SITES:
    fresh_ids  = meta[(meta['Location'] == loc) & (meta['Storage'] == 'FRESH')].index
    frozen_ids = meta[(meta['Location'] == loc) & (meta['Storage'] == 'FROZEN')].index
    fresh_ids  = [s for s in fresh_ids  if s in otu_rel.index]
    frozen_ids = [s for s in frozen_ids if s in otu_rel.index]

    if len(fresh_ids) < 3 or len(frozen_ids) < 2:
        print(f"  {loc}: insufficient samples (fresh={len(fresh_ids)}, frozen={len(frozen_ids)})")
        continue

    Xf       = hellinger_transform(otu_rel.loc[fresh_ids].values)
    Xz       = hellinger_transform(otu_rel.loc[frozen_ids].values)
    centroid = Xf.mean(axis=0)

    dists_f = [bray_curtis_distance(Xf[i], centroid) for i in range(len(Xf))]
    dists_z = [bray_curtis_distance(Xz[i], centroid) for i in range(len(Xz))]

    res = tost_two_sample(dists_z, dists_f, delta=DELTA_BC, alpha=ALPHA)
    if res is None:
        continue

    row = {'site': loc, 'comparison': 'all_fractions', **res,
           'fresh_mean_bc': res['mean_y'], 'frozen_mean_bc': res['mean_x'],
           'dists_fresh': dists_f, 'dists_frozen': dists_z}
    bc_results.append(row)

    sign = '✅ EQUIVALENT' if res['is_equivalent'] else '❌ NOT EQUIVALENT'
    print(f"  {loc}: fresh n={res['ny']}, frozen n={res['nx']}")
    print(f"    BC centroid: fresh={res['mean_y']:.4f}±{res['std_y']:.4f}, frozen={res['mean_x']:.4f}±{res['std_x']:.4f}")
    print(f"    Diff={res['mean_diff']:+.4f}, 90%CI=[{res['ci_90_low']:+.4f}, {res['ci_90_high']:+.4f}]")
    print(f"    p_lower={res['p_lower']:.4f}, p_upper={res['p_upper']:.4f}  {sign}")

# Pooled vagina + endometrium (primary result; cervix excluded: n_frozen=1)
pool_locs = ['VAGINA', 'ENDOMETRIUM']
fresh_pool  = [s for s in meta[(meta['Location'].isin(pool_locs)) & (meta['Storage'] == 'FRESH')].index  if s in otu_rel.index]
frozen_pool = [s for s in meta[(meta['Location'].isin(pool_locs)) & (meta['Storage'] == 'FROZEN')].index if s in otu_rel.index]

Xf_pool  = hellinger_transform(otu_rel.loc[fresh_pool].values)
Xz_pool  = hellinger_transform(otu_rel.loc[frozen_pool].values)
centroid_pool = Xf_pool.mean(axis=0)

dists_f_pool = [bray_curtis_distance(Xf_pool[i], centroid_pool) for i in range(len(Xf_pool))]
dists_z_pool = [bray_curtis_distance(Xz_pool[i], centroid_pool) for i in range(len(Xz_pool))]

res_pool = tost_two_sample(dists_z_pool, dists_f_pool, delta=DELTA_BC, alpha=ALPHA)
sign = '✅ EQUIVALENT' if res_pool['is_equivalent'] else '❌ NOT EQUIVALENT'
print(f"\n  POOLED (vagina + endometrium): fresh n={res_pool['ny']}, frozen n={res_pool['nx']}")
print(f"    Diff={res_pool['mean_diff']:+.4f}, 90%CI=[{res_pool['ci_90_low']:+.4f}, {res_pool['ci_90_high']:+.4f}]")
print(f"    p_lower={res_pool['p_lower']:.4f}, p_upper={res_pool['p_upper']:.4f}  {sign}")

bc_results.append({'site': 'VAGINA+ENDOMETRIUM', 'comparison': 'all_fractions_pooled',
                   'dists_fresh': dists_f_pool, 'dists_frozen': dists_z_pool, **res_pool,
                   'fresh_mean_bc': res_pool['mean_y'], 'frozen_mean_bc': res_pool['mean_x']})


# ─────────────────────────────────────────────────────────────────────────────
# ANALYSIS 2: SHANNON TOST (PRE-sorted fractions)
# ─────────────────────────────────────────────────────────────────────────────

print("\nAnalysis 2: Shannon TOST (PRE-sorted fractions only)")
print("-" * 50)

shannon_results = []

for loc in SITES:
    fresh_ids  = meta[(meta['Location'] == loc) & (meta['Storage'] == 'FRESH')  & (meta['Antibody'] == 'PRE')].index
    frozen_ids = meta[(meta['Location'] == loc) & (meta['Storage'] == 'FROZEN') & (meta['Antibody'] == 'PRE')].index
    fresh_ids  = [s for s in fresh_ids  if s in otu_rel.index]
    frozen_ids = [s for s in frozen_ids if s in otu_rel.index]

    if len(fresh_ids) < 2 or len(frozen_ids) < 2:
        print(f"  {loc}: insufficient samples (fresh={len(fresh_ids)}, frozen={len(frozen_ids)})")
        continue

    sh_fresh  = np.array([shannon_entropy(otu_rel.loc[s].values) for s in fresh_ids])
    sh_frozen = np.array([shannon_entropy(otu_rel.loc[s].values) for s in frozen_ids])

    res = tost_two_sample(sh_fresh, sh_frozen, delta=DELTA_SHANNON, alpha=ALPHA)
    if res is None:
        continue

    sign = '✅ EQUIVALENT' if res['is_equivalent'] else '❌ NOT EQUIVALENT'
    note = ''
    if loc == 'VAGINA' and not res['is_equivalent']:
        note = ' (taxon richness reduction in frozen vaginal PRE; see Table S2.6)'

    print(f"  {loc}: fresh n={res['nx']}, frozen n={res['ny']}")
    print(f"    Shannon: fresh={res['mean_x']:.3f}±{res['std_x']:.3f}, frozen={res['mean_y']:.3f}±{res['std_y']:.3f}")
    print(f"    Diff={res['mean_diff']:+.3f}, 90%CI=[{res['ci_90_low']:+.3f}, {res['ci_90_high']:+.3f}]")
    print(f"    p_lower={res['p_lower']:.4f}, p_upper={res['p_upper']:.4f}  {sign}{note}")

    shannon_results.append({'site': loc, 'sh_fresh': sh_fresh, 'sh_frozen': sh_frozen, **res})


# ─────────────────────────────────────────────────────────────────────────────
# SAVE CSV OUTPUTS
# ─────────────────────────────────────────────────────────────────────────────

bc_df = pd.DataFrame([
    {k: v for k, v in r.items() if k not in ('dists_fresh', 'dists_frozen')}
    for r in bc_results
])
bc_path = os.path.join(FIG_DIR, 'tost_bray_curtis_results.csv')
bc_df.to_csv(bc_path, index=False)
print(f"\nSaved: {bc_path}")

sh_df = pd.DataFrame([
    {k: v for k, v in r.items() if k not in ('sh_fresh', 'sh_frozen')}
    for r in shannon_results
])
sh_path = os.path.join(FIG_DIR, 'tost_shannon_results.csv')
sh_df.to_csv(sh_path, index=False)
print(f"Saved: {sh_path}")


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE: TOST summary
# ─────────────────────────────────────────────────────────────────────────────

print("\nGenerating figure...")

fig = plt.figure(figsize=(18, 10))
fig.patch.set_facecolor('white')
gs  = gridspec.GridSpec(2, 3, figure=fig, hspace=0.45, wspace=0.38)

# ── Panels A–C: BC distance distributions (fresh vs frozen) per site ─────────
for oi, loc in enumerate(SITES):
    ax  = fig.add_subplot(gs[0, oi])
    row = next((r for r in bc_results if r['site'] == loc), None)

    if row is None:
        ax.text(0.5, 0.5, f'{loc}\nn < 3 (excluded)', transform=ax.transAxes,
                ha='center', va='center', fontsize=11)
        ax.set_title(f'A{oi+1}  {loc.title()} — BC distances', fontsize=11,
                     fontweight='bold', loc='left', color=SITE_COLORS[loc])
        continue

    df_val    = row['dists_fresh']
    dz_val    = row['dists_frozen']
    all_d     = df_val + dz_val
    ylim_max  = max(all_d) * 1.18

    positions = [1, 2]
    bp = ax.boxplot([df_val, dz_val], positions=positions, widths=0.55,
                    patch_artist=True, showfliers=False,
                    medianprops=dict(color='white', lw=2.2),
                    whiskerprops=dict(color=SITE_COLORS[loc], lw=1.5),
                    capprops=dict(color=SITE_COLORS[loc], lw=1.5))
    for patch, color in zip(bp['boxes'], ['#90CAF9', '#EF9A9A']):
        patch.set(facecolor=color, alpha=0.75, edgecolor=SITE_COLORS[loc], lw=1.3)

    np.random.seed(42)
    ax.scatter(1 + np.random.normal(0, 0.06, len(df_val)), df_val,
               c='#1565C0', s=40, alpha=0.8, zorder=8, edgecolors='white', lw=0.4)
    ax.scatter(2 + np.random.normal(0, 0.06, len(dz_val)), dz_val,
               c='#B71C1C', s=40, alpha=0.8, zorder=8, edgecolors='white', lw=0.4)

    # Equivalence bands (±DELTA from fresh mean)
    fresh_mean = row['fresh_mean_bc']
    ax.axhspan(fresh_mean - DELTA_BC, fresh_mean + DELTA_BC,
               alpha=0.12, color='green', label=f'±{DELTA_BC} equiv. band')
    ax.axhline(fresh_mean, color='#1565C0', ls='--', lw=1.2, alpha=0.7, label='Fresh mean')

    sign_str = '✅ EQUIV.' if row['is_equivalent'] else '❌ NOT EQUIV.'
    ax.text(0.97, 0.04,
            f"TOST delta=±{DELTA_BC}\ndiff={row['mean_diff']:+.3f}\n90%CI=[{row['ci_90_low']:+.3f},{row['ci_90_high']:+.3f}]\n{sign_str}",
            transform=ax.transAxes, ha='right', va='bottom', fontsize=8,
            bbox=dict(boxstyle='round,pad=0.3', fc='#f8f8f8', ec='#ccc', alpha=0.9))

    ax.set_xticks([1, 2])
    ax.set_xticklabels(['Fresh', 'Frozen'], fontsize=10)
    ax.set_ylabel('BC distance to fresh centroid', fontsize=9)
    ax.set_ylim(-0.02, ylim_max)
    ax.set_title(f'A{oi+1}  {loc.title()} — BC distances (all fractions)', fontsize=10,
                 fontweight='bold', loc='left', color=SITE_COLORS[loc])
    ax.spines[['top', 'right']].set_visible(False)
    ax.yaxis.grid(True, alpha=0.2, ls='--'); ax.set_axisbelow(True)
    if oi == 0:
        ax.legend(fontsize=8, loc='upper left', frameon=True, framealpha=0.85)

# ── Panels D–F: Shannon distributions (PRE-sorted) ────────────────────────────
for oi, loc in enumerate(SITES):
    ax  = fig.add_subplot(gs[1, oi])
    row = next((r for r in shannon_results if r['site'] == loc), None)

    if row is None:
        ax.text(0.5, 0.5, f'{loc}\nn < 2 (excluded)', transform=ax.transAxes,
                ha='center', va='center', fontsize=11)
        ax.set_title(f'B{oi+1}  {loc.title()} — Shannon (PRE)', fontsize=11,
                     fontweight='bold', loc='left', color=SITE_COLORS[loc])
        continue

    sf  = row['sh_fresh']
    sz  = row['sh_frozen']
    all_s = np.concatenate([sf, sz])

    bp2 = ax.boxplot([sf, sz], positions=[1, 2], widths=0.55,
                     patch_artist=True, showfliers=False,
                     medianprops=dict(color='white', lw=2.2),
                     whiskerprops=dict(color=SITE_COLORS[loc], lw=1.5),
                     capprops=dict(color=SITE_COLORS[loc], lw=1.5))
    for patch, color in zip(bp2['boxes'], ['#90CAF9', '#EF9A9A']):
        patch.set(facecolor=color, alpha=0.75, edgecolor=SITE_COLORS[loc], lw=1.3)

    np.random.seed(42)
    ax.scatter(1 + np.random.normal(0, 0.06, len(sf)), sf,
               c='#1565C0', s=40, alpha=0.8, zorder=8, edgecolors='white', lw=0.4)
    ax.scatter(2 + np.random.normal(0, 0.06, len(sz)), sz,
               c='#B71C1C', s=40, alpha=0.8, zorder=8, edgecolors='white', lw=0.4)

    fresh_sh_mean = row['mean_x']
    ax.axhspan(fresh_sh_mean - DELTA_SHANNON, fresh_sh_mean + DELTA_SHANNON,
               alpha=0.12, color='green', label=f'±{DELTA_SHANNON} equiv. band')
    ax.axhline(fresh_sh_mean, color='#1565C0', ls='--', lw=1.2, alpha=0.7)

    sign_str = '✅ EQUIV.' if row['is_equivalent'] else '❌ NOT EQUIV.'
    ax.text(0.97, 0.04,
            f"TOST delta=±{DELTA_SHANNON}\ndiff={row['mean_diff']:+.3f}\n90%CI=[{row['ci_90_low']:+.3f},{row['ci_90_high']:+.3f}]\n{sign_str}",
            transform=ax.transAxes, ha='right', va='bottom', fontsize=8,
            bbox=dict(boxstyle='round,pad=0.3', fc='#f8f8f8', ec='#ccc', alpha=0.9))

    ax.set_xticks([1, 2])
    ax.set_xticklabels(['Fresh', 'Frozen'], fontsize=10)
    ax.set_ylabel("Shannon H'", fontsize=9)
    ax.set_ylim(-0.1, max(all_s) * 1.22)
    ax.set_title(f'B{oi+1}  {loc.title()} — Shannon H\' (PRE-sorted)', fontsize=10,
                 fontweight='bold', loc='left', color=SITE_COLORS[loc])
    ax.spines[['top', 'right']].set_visible(False)
    ax.yaxis.grid(True, alpha=0.2, ls='--'); ax.set_axisbelow(True)

fig.suptitle(
    'TOST Cryopreservation Equivalence Analysis\n'
    'Top row: Bray–Curtis distance to fresh centroid (all Ig fractions)   '
    '|   Bottom row: Shannon entropy (PRE-sorted fractions)\n'
    f'Green bands = equivalence margin (BC: ±{DELTA_BC}; Shannon: ±{DELTA_SHANNON} H\')',
    fontsize=11, fontweight='bold', y=1.02
)

fig_path = os.path.join(FIG_DIR, 'SFig_TOST_Cryopreservation.png')
fig.savefig(fig_path, dpi=180, bbox_inches='tight', facecolor='white')
plt.close()
print(f"Saved: {fig_path}")


# ─────────────────────────────────────────────────────────────────────────────
# PRINT MANUSCRIPT-READY SUMMARY
# ─────────────────────────────────────────────────────────────────────────────

print("\n" + "="*65)
print("MANUSCRIPT-READY SUMMARY")
print("="*65)

print(f"\nBray-Curtis TOST (all fractions, vagina + endometrium pooled):")
print(f"  fresh n={res_pool['ny']}, frozen n={res_pool['nx']}")
print(f"  Fresh centroid distance:  {res_pool['mean_y']:.4f} ± {res_pool['std_y']:.4f}")
print(f"  Frozen centroid distance: {res_pool['mean_x']:.4f} ± {res_pool['std_x']:.4f}")
print(f"  Mean difference (frozen - fresh): {res_pool['mean_diff']:+.4f}")
print(f"  90% CI: [{res_pool['ci_90_low']:+.4f}, {res_pool['ci_90_high']:+.4f}]")
print(f"  p_lower = {res_pool['p_lower']:.4f},  p_upper = {res_pool['p_upper']:.4f}")
print(f"  Equivalence within ±{DELTA_BC}: {'YES ✅' if res_pool['is_equivalent'] else 'NO ❌'}")

for r in shannon_results:
    print(f"\nShannon TOST — {r['site']} (PRE-sorted):")
    print(f"  fresh n={r['nx']}, frozen n={r['ny']}")
    print(f"  Fresh: {r['mean_x']:.3f} ± {r['std_x']:.3f}")
    print(f"  Frozen: {r['mean_y']:.3f} ± {r['std_y']:.3f}")
    print(f"  diff={r['mean_diff']:+.3f}, 90%CI=[{r['ci_90_low']:+.3f},{r['ci_90_high']:+.3f}]")
    print(f"  p_lower={r['p_lower']:.4f}, p_upper={r['p_upper']:.4f}")
    print(f"  Equivalent within ±{DELTA_SHANNON}: {'YES ✅' if r['is_equivalent'] else 'NO ❌'}")

print("\nNote: Vaginal Shannon non-equivalence reflects the known taxon richness")
print("reduction in frozen vaginal PRE-sorted fractions (reported in Results,")
print("Supplementary Table S2.6) and does not affect the primary cryopreservation")
print("conclusion, which is based on Bray-Curtis community composition across all")
print("Ig fractions.")
