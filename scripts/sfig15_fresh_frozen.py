import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Ellipse
from scipy.stats import mannwhitneyu
import warnings; warnings.filterwarnings('ignore')

from shared_utils import (load_data, shannon, hellinger, bray_curtis, pcoa,
                           permanova, pct_label,
                           LOC_C, STOR_C, LOC_ORDER, AB_ORDER, AB_LABELS, FIG_DIR)

otu_rel, meta, tax = load_data()
DATA  = otu_rel.T
taxa  = otu_rel.index.tolist()

def plot_ellipse(ax, pts, color, n_std=1.96, alpha=0.15):
    if len(pts) < 3: return
    mean = pts.mean(axis=0); cov = np.cov(pts.T)
    if cov.ndim < 2 or np.any(np.isnan(cov)): return
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]; vals = vals[order]; vecs = vecs[:, order]
    angle = np.degrees(np.arctan2(*vecs[:, 0][::-1]))
    w, h = 2 * n_std * np.sqrt(np.abs(vals))
    ax.add_patch(Ellipse(xy=mean, width=w, height=h, angle=angle,
                         facecolor=color, alpha=alpha, edgecolor=color, lw=1.5))

# ── Diversity table ───────────────────────────────────────────────────────────
div = pd.DataFrame([{'SampleID': s, 'Shannon': shannon(DATA.loc[s].values)}
                    for s in DATA.index]).set_index('SampleID').join(meta)

# ─────────────────────────────────────────────────────────────────────────────
# MAIN FIGURE
# ─────────────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(22, 18))
fig.patch.set_facecolor('white')

gs = fig.add_gridspec(2, 3, hspace=0.15, wspace=0.15)

# ── Panels A: PCoA per organ, coloured by storage ────────────────────────────
for oi, loc in enumerate(LOC_ORDER):
    ax   = fig.add_subplot(gs[0, oi])
    sids = meta[meta['Location'] == loc].index
    sids = [s for s in sids if s in DATA.index]
    if len(sids) < 4: continue

    Xi   = otu_rel[sids].T.values
    Xh   = hellinger(Xi)
    Di   = bray_curtis(Xh)
    ci_c, vi = pcoa(Di, n_axes=2)
    stor = meta.loc[sids, 'Storage'].values

    for s in ['FRESH', 'FROZEN']:
        idx = np.where(stor == s)[0]
        pts = ci_c[idx, :2]
        col = STOR_C[s]
        ax.scatter(pts[:, 0], pts[:, 1], c=col, s=65, alpha=0.85,
                   edgecolors='white', lw=0.5, label=s.title(), zorder=5)
        plot_ellipse(ax, pts, col)

    pr = permanova(Di, stor, n_perm=499)
    ax.text(0.97, 0.03,
            f"PERMANOVA\nR²={pr['R2']:.3f}, p={pr['p']:.3f}",
            transform=ax.transAxes, ha='right', va='bottom', fontsize=9,
            style='italic',
            bbox=dict(boxstyle='round,pad=0.3', fc='#f8f8f8', ec='#ccc', alpha=0.85))
    ax.set_title(f'A   {loc.title()} — PCoA by storage', fontsize=11,
                 fontweight='bold', color=LOC_C[loc], loc='left', pad=6)
    ax.set_xlabel(f'PCoA1 ({vi[0]*100:.1f}%)', fontsize=9)
    ax.set_ylabel(f'PCoA2 ({vi[1]*100:.1f}%)', fontsize=9)
    ax.spines[['top', 'right']].set_visible(False)
    ax.axhline(0, color='#ddd', lw=0.8, linestyle='--')
    ax.axvline(0, color='#ddd', lw=0.8, linestyle='--')
    _n_fresh_a  = sum(1 for s in sids if meta.loc[s, 'Storage'] == 'FRESH')
    _n_frozen_a = sum(1 for s in sids if meta.loc[s, 'Storage'] == 'FROZEN')
    ax.legend(handles=[
        mpatches.Patch(color=STOR_C['FRESH'],  label=f"Fresh (n={_n_fresh_a})"),
        mpatches.Patch(color=STOR_C['FROZEN'], label=f"Frozen (n={_n_frozen_a})"),
    ], fontsize=9, frameon=True, loc='lower center', ncol=2, bbox_to_anchor=(0.5, -0.12))

# ── Panels B: Shannon Fresh vs Frozen per organ ───────────────────────────────
for oi, loc in enumerate(LOC_ORDER):
    ax   = fig.add_subplot(gs[1, oi])
    sub  = div[div['Location'] == loc]

    xtick_pos, xtick_lab = [], []
    _all_sh = sub['Shannon'].dropna()
    _y_span = max(_all_sh.max() - _all_sh.min(), 0.1) if len(_all_sh) > 0 else 1.0
    _y_pad  = 0.08 * _y_span
    _y_lift = 0.03 * _y_span
    _y_text = 0.02 * _y_span

    for ai, ab in enumerate(AB_ORDER):
        for si, stor in enumerate(['FRESH', 'FROZEN']):
            grp = sub[(sub['Antibody'] == ab) & (sub['Storage'] == stor)]['Shannon'].dropna()
            if len(grp) == 0: continue
            col = STOR_C[stor]
            xp  = ai * 2.7 + si * 0.9
            ax.boxplot(grp, positions=[xp], widths=0.75, patch_artist=True,
                       showfliers=False,
                       medianprops={'color': 'white', 'lw': 2},
                       whiskerprops={'color': col, 'lw': 1.5},
                       capprops={'color': col, 'lw': 1.5},
                       boxprops={'facecolor': col, 'alpha': 0.50,
                                  'edgecolor': col, 'lw': 1.2})
            np.random.seed(ai * 10 + si)
            ax.scatter(xp + np.random.normal(0, 0.08, len(grp)), grp,
                       c=col, s=35, alpha=0.8, zorder=8, edgecolors='white', lw=0.3)
        xtick_pos.append(ai * 2.7 + 0.45)
        xtick_lab.append(AB_LABELS[ab])

        # MWU
        h = sub[(sub['Antibody'] == ab) & (sub['Storage'] == 'FRESH')]['Shannon'].dropna()
        f_v = sub[(sub['Antibody'] == ab) & (sub['Storage'] == 'FROZEN')]['Shannon'].dropna()
        if len(h) > 1 and len(f_v) > 1:
            try:
                _, p = mannwhitneyu(h, f_v, alternative='two-sided')
                sig = pct_label(p)
                x1  = ai * 2.7
                x2  = ai * 2.7 + 0.9
                _ym = max(h.max(), f_v.max()) + _y_pad
                ax.plot([x1, x1, x2, x2],
                        [_ym, _ym + _y_lift, _ym + _y_lift, _ym],
                        color='black', lw=0.9, zorder=11)
                ax.text((x1 + x2) / 2, _ym + _y_lift + _y_text, sig,
                        ha='center', va='bottom', fontsize=8,
                        fontweight='bold', color='black')
            except: pass

    ax.set_xticks(xtick_pos)
    ax.set_xticklabels(xtick_lab, fontsize=9, rotation=0, ha='center')
    ax.set_ylabel("Shannon H'", fontsize=10)
    ax.set_title(f'B   {loc.title()} — Shannon by storage', fontsize=11,
                 fontweight='bold', color=LOC_C[loc], loc='left', pad=6)
    ax.spines[['top', 'right']].set_visible(False)
    ax.yaxis.grid(True, alpha=0.22, linestyle='--')
    ax.set_axisbelow(True)
    _n_fresh_b  = len(sub[sub['Storage'] == 'FRESH'])
    _n_frozen_b = len(sub[sub['Storage'] == 'FROZEN'])
    ax.legend(handles=[
        mpatches.Patch(color=STOR_C['FRESH'],  label=f"Fresh (n={_n_fresh_b})"),
        mpatches.Patch(color=STOR_C['FROZEN'], label=f"Frozen (n={_n_frozen_b})"),
    ], fontsize=9, frameon=True, loc='lower center',
       bbox_to_anchor=(0.5, -0.1), ncol=2)

plt.tight_layout()
fig.savefig(os.path.join(FIG_DIR, 'SupplementaryFigure15.png'),
            dpi=600, bbox_inches='tight', facecolor='white')
plt.close()
print("Saved SupplementaryFigure15.png")