import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Ellipse
import warnings; warnings.filterwarnings('ignore')

from shared_utils import (load_data, hellinger, bray_curtis, pcoa,
                           permanova,
                           LOC_C, AB_C,
                           LOC_ORDER, AB_ORDER, AB_LABELS, FIG_DIR)

otu_rel, meta, tax = load_data()


# ── Confidence ellipse helper ─────────────────────────────────────────────────
def plot_ellipse(ax, pts, color, alpha=0.18, n_std=1.96):
    if len(pts) < 3:
        return
    mean  = pts.mean(axis=0)
    cov   = np.cov(pts.T)
    if cov.ndim < 2 or np.any(np.isnan(cov)):
        return
    vals, vecs = np.linalg.eigh(cov)
    order  = vals.argsort()[::-1]
    vals   = vals[order]; vecs = vecs[:, order]
    angle  = np.degrees(np.arctan2(*vecs[:, 0][::-1]))
    w, h   = 2 * n_std * np.sqrt(np.abs(vals))
    ell    = Ellipse(xy=mean, width=w, height=h, angle=angle,
                     facecolor=color, alpha=alpha, edgecolor=color, linewidth=1.5)
    ax.add_patch(ell)


# ─────────────────────────────────────────────────────────────────────────────
# SUPPLEMENTARY: PCoA per location, coloured by Antibody fraction
# ─────────────────────────────────────────────────────────────────────────────
fig2, axes2 = plt.subplots(1, 3, figsize=(20, 7))
fig2.patch.set_facecolor('white')

for oi, loc in enumerate(LOC_ORDER):
    ax   = axes2[oi]
    sids = meta[meta['Location'] == loc].index.tolist()
    sids = [s for s in sids if s in otu_rel.columns]
    if len(sids) < 4:
        ax.text(0.5, 0.5, 'Insufficient data', transform=ax.transAxes,
                ha='center', va='center')
        continue

    Xi   = otu_rel[sids].T.values
    Xh   = hellinger(Xi)
    Di   = bray_curtis(Xh)
    ci_c, vi = pcoa(Di, n_axes=2)

    ab_labels = meta.loc[sids, 'Antibody'].values
    for ab in AB_ORDER:
        idx  = np.where(ab_labels == ab)[0]
        pts  = ci_c[idx, :2]
        col  = AB_C[ab]
        ax.scatter(pts[:, 0], pts[:, 1], c=col, s=70, alpha=0.85,
                   edgecolors='white', lw=0.5, label=AB_LABELS[ab], zorder=5)
        plot_ellipse(ax, pts, col)

    # PERMANOVA
    pr = permanova(Di, ab_labels, n_perm=499)
    ax.text(0.97, 0.03,
            f"PERMANOVA\nR²={pr['R2']:.3f}, p={pr['p']:.3f}",
            transform=ax.transAxes, ha='right', va='bottom', fontsize=9,
            style='italic', color='#333',
            bbox=dict(boxstyle='round,pad=0.3', fc='#f8f8f8', ec='#ccc', alpha=0.85))

    ax.set_title(loc.title(), fontsize=13, fontweight='bold',
                 color=LOC_C[loc], pad=8)
    ax.set_xlabel(f'PCoA1 ({vi[0]*100:.1f}%)', fontsize=10)
    ax.set_ylabel(f'PCoA2 ({vi[1]*100:.1f}%)' if oi == 0 else '', fontsize=10)
    ax.spines[['top', 'right']].set_visible(False)
    ax.axhline(0, color='#ccc', lw=0.8, linestyle='--')
    ax.axvline(0, color='#ccc', lw=0.8, linestyle='--')

handles = [mpatches.Patch(color=AB_C[ab], label=AB_LABELS[ab]) for ab in AB_ORDER]
fig2.legend(handles=handles, loc='lower center', ncol=4, fontsize=11,
            frameon=True, bbox_to_anchor=(0.5, -0.04))

plt.tight_layout()
fig2.savefig(os.path.join(FIG_DIR, 'SupplementaryFigure3.png'),
             dpi=600, bbox_inches='tight', facecolor='white')
plt.close()
print("Saved SupplementaryFigure3.png")
