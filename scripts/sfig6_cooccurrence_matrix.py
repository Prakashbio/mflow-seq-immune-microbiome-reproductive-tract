import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
import warnings; warnings.filterwarnings('ignore')

from shared_utils import (load_data, LOC_C, LOC_ORDER, FIG_DIR)

otu_rel, meta, tax = load_data()
DATA  = otu_rel.T
taxa  = otu_rel.index.tolist()

THRESH_R = 0.4

# ─────────────────────────────────────────────────────────────────────────────
# SUPPLEMENTARY: Full co-occurrence matrix heatmap
# ─────────────────────────────────────────────────────────────────────────────
fig2, axes2 = plt.subplots(1, 3, figsize=(24, 8))
fig2.patch.set_facecolor('white')

for oi, loc in enumerate(LOC_ORDER):
    ax   = axes2[oi]
    sids = meta[(meta['Location'] == loc) & (meta['Antibody'] == 'PRE')].index
    sids = [s for s in sids if s in DATA.index]
    if not sids: continue

    sub = DATA.loc[sids]
    n   = len(taxa)
    cm  = np.zeros((n, n))
    pm  = np.ones((n, n))
    for i in range(n):
        for j in range(i, n):
            xi = sub.iloc[:, i].values
            xj = sub.iloc[:, j].values
            if np.std(xi) > 0 and np.std(xj) > 0:
                r, p = spearmanr(xi, xj)
                cm[i, j] = cm[j, i] = r
                pm[i, j] = pm[j, i] = p
            if i == j: cm[i, j] = 1.0

    im = ax.imshow(cm, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto',
                   interpolation='nearest')
    # Mark significant cells
    sig_i, sig_j = np.where((pm < 0.05) & (np.abs(cm) > THRESH_R) & (~np.eye(n, dtype=bool)))
    for si, sj in zip(sig_i, sig_j):
        ax.add_patch(plt.Rectangle((sj - 0.5, si - 0.5), 1, 1,
                                    fill=False, edgecolor='black', lw=0.8))

    ax.set_xticks(range(n))
    ax.set_xticklabels(taxa, rotation=90, ha='right', fontsize=6.5, style='italic')
    ax.set_yticks(range(n))
    ax.set_yticklabels(taxa, fontsize=6.5, style='italic')
    ax.set_title(loc.title(), fontsize=13, fontweight='bold', color=LOC_C[loc], pad=8)
    plt.colorbar(im, ax=ax, label='Spearman r', fraction=0.04, pad=0.02)

fig2.text(0.5, -0.02, 'Black borders = significant (|r|≥0.4, p<0.05)',
          ha='center', fontsize=10, style='italic')

plt.tight_layout()
fig2.savefig(os.path.join(FIG_DIR, 'SupplementaryFigure6.png'),
             dpi=600, bbox_inches='tight', facecolor='white')
plt.close()
print("Saved SupplementaryFigure6.png")
