
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats import spearmanr
import warnings; warnings.filterwarnings('ignore')

from shared_utils import (load_data, LOC_C, LOC_ORDER, FIG_DIR)

otu_rel, meta, tax = load_data()
DATA  = otu_rel.T
taxa  = otu_rel.index.tolist()

# ─────────────────────────────────────────────────────────────────────────────
# SUPPLEMENTARY: Per-taxon Fresh vs Frozen scatter correlation
# ─────────────────────────────────────────────────────────────────────────────
fig2, axes2 = plt.subplots(1, 3, figsize=(20, 7))
fig2.patch.set_facecolor('white')

for oi, loc in enumerate(LOC_ORDER):
    ax = axes2[oi]

    fresh_sids  = meta[(meta['Location'] == loc) & (meta['Storage'] == 'FRESH')].index
    frozen_sids = meta[(meta['Location'] == loc) & (meta['Storage'] == 'FROZEN')].index
    fresh_sids  = [s for s in fresh_sids  if s in DATA.index]
    frozen_sids = [s for s in frozen_sids if s in DATA.index]

    if not fresh_sids or not frozen_sids:
        ax.text(0.5, 0.5, 'Insufficient data', transform=ax.transAxes, ha='center')
        continue

    mu_fresh  = DATA.loc[fresh_sids].mean(axis=0).values * 100
    mu_frozen = DATA.loc[frozen_sids].mean(axis=0).values * 100

    r, p = spearmanr(mu_fresh, mu_frozen)

    ax.scatter(mu_fresh, mu_frozen, s=80, alpha=0.80,
               c=LOC_C[loc], edgecolors='white', lw=0.5, zorder=5)

    # diagonal
    maxv = max(mu_fresh.max(), mu_frozen.max()) * 1.05
    ax.plot([0, maxv], [0, maxv], 'k--', lw=1.2, alpha=0.35)

    for i, t in enumerate(taxa):
        if mu_fresh[i] + mu_frozen[i] > 3:
            ax.annotate(t.split(' ')[1], xy=(mu_fresh[i], mu_frozen[i]),
                        xytext=(4, 3), textcoords='offset points',
                        fontsize=7.5, style='italic', color='#333')

    _p_lbl = ('p < 0.001' if p < 0.001 else 'p < 0.01' if p < 0.01
              else 'p < 0.05' if p < 0.05 else 'p ≥ 0.05')
    ax.text(0.05, 0.92, f'Spearman r = {r:.3f}\n{_p_lbl}',
            transform=ax.transAxes, fontsize=10, color='#333',
            bbox=dict(boxstyle='round,pad=0.3', fc='#f8f8f8', ec='#ccc', alpha=0.85))

    ax.set_xlabel('Mean relative abundance — Fresh (%)', fontsize=10)
    ax.set_ylabel('Mean relative abundance — Frozen (%)', fontsize=10)
    ax.set_title(loc.title(), fontsize=13, fontweight='bold', color=LOC_C[loc], pad=8)
    ax.spines[['top', 'right']].set_visible(False)
    ax.grid(True, alpha=0.18, linestyle='--')

plt.tight_layout()
fig2.savefig(os.path.join(FIG_DIR, 'SupplementaryFigure7.png'),
             dpi=600, bbox_inches='tight', facecolor='white')
plt.close()
print("Saved SupplementaryFigure7.png")
