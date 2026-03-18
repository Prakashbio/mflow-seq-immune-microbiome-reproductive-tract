import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings; warnings.filterwarnings('ignore')

from shared_utils import (load_data, LOC_C, LOC_ORDER, FIG_DIR)

otu_rel, meta, tax = load_data()
DATA  = otu_rel.T
taxa  = otu_rel.index.tolist()

# ─────────────────────────────────────────────────────────────────────────────
# SUPPLEMENTARY: Taxon accumulation curve
# ─────────────────────────────────────────────────────────────────────────────
fig2, axes2 = plt.subplots(1, 3, figsize=(20, 7))
fig2.patch.set_facecolor('white')

N_PERM = 50

for oi, loc in enumerate(LOC_ORDER):
    ax   = axes2[oi]
    sids = meta[meta['Location'] == loc].index
    sids = [s for s in sids if s in DATA.index]
    if not sids: continue

    sub    = DATA.loc[sids]
    n_samp = len(sids)
    accum  = np.zeros((N_PERM, n_samp))

    np.random.seed(42)
    for pi in range(N_PERM):
        order = np.random.permutation(n_samp)
        seen  = set()
        for si, idx in enumerate(order):
            seen.update(np.where(sub.iloc[idx].values > 0)[0])
            accum[pi, si] = len(seen)

    mean_acc = accum.mean(axis=0)
    lo_acc   = np.percentile(accum, 10, axis=0)
    hi_acc   = np.percentile(accum, 90, axis=0)
    x        = np.arange(1, n_samp + 1)

    ax.fill_between(x, lo_acc, hi_acc, alpha=0.25, color=LOC_C[loc])
    ax.plot(x, mean_acc, color=LOC_C[loc], lw=2.5, label='Mean accumulation')
    ax.axhline(len(taxa), color='#888', lw=1.5, linestyle=':', alpha=0.8,
               label=f'Total taxa ({len(taxa)})')

    ax.set_xlabel('Number of samples', fontsize=10)
    ax.set_ylabel('Cumulative taxa observed', fontsize=10)
    ax.set_title(loc.title(), fontsize=13, fontweight='bold', color=LOC_C[loc], pad=8)
    ax.legend(fontsize=9, frameon=True)
    ax.spines[['top', 'right']].set_visible(False)
    ax.grid(True, alpha=0.18, linestyle='--')
    ax.set_xlim(1, n_samp)
    ax.set_ylim(0, len(taxa) + 3)

plt.tight_layout()
fig2.savefig(os.path.join(FIG_DIR, 'SupplementaryFigure8.png'),
             dpi=600, bbox_inches='tight', facecolor='white')
plt.close()
print("Saved SupplementaryFigure8.png")
