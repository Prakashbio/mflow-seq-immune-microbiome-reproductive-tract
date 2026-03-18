import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import warnings; warnings.filterwarnings('ignore')
from shared_utils import (load_data, LOC_C, LOC_ORDER, AB_LABELS, FIG_DIR)

otu_rel, meta, tax = load_data()
DATA   = otu_rel.T
taxa   = otu_rel.index.tolist()
PSEUDO = 1e-5

# ── Compute mean relative abundance per (location, fraction, taxon) ───────────
def group_mean(loc, ab):
    sids = meta[(meta['Location'] == loc) & (meta['Antibody'] == ab)].index
    sids = [s for s in sids if s in DATA.index]
    if not sids:
        return pd.Series(0.0, index=taxa)
    return DATA.loc[sids].mean(axis=0)


# ─────────────────────────────────────────────────────────────────────────────
# SUPPLEMENTARY: IgA vs Pre-sorted and IgM vs Pre-sorted log2FC
# ─────────────────────────────────────────────────────────────────────────────
fig2, axes2 = plt.subplots(2, 3, figsize=(22, 16))
fig2.patch.set_facecolor('white')

for oi, loc in enumerate(LOC_ORDER):
    mu_pre = group_mean(loc, 'PRE')
    for ri, (ab, ab_label, col) in enumerate([('A', 'IgA', '#E66101'), ('M', 'IgM', '#5E3C99')]):
        ax  = axes2[ri, oi]
        mu  = group_mean(loc, ab)
        lfc = np.log2((mu.values + PSEUDO) / (mu_pre.values + PSEUDO))
        df  = pd.DataFrame({'Taxon': taxa, 'log2FC': lfc}).sort_values('log2FC')
        y   = np.arange(len(df))
        colors_bar = [col if v > 0 else '#AAAAAA' for v in df['log2FC']]
        ax.barh(y, df['log2FC'], color=colors_bar, alpha=0.72, height=0.65, edgecolor='white')
        ax.axvline(0, color='#333', lw=1.5, linestyle='--', alpha=0.7)
        ax.set_yticks(y)
        ax.set_yticklabels(df['Taxon'], fontsize=7.5, style='italic')
        ax.set_xlabel(f'log₂({ab_label} / Pre-sorted)', fontsize=9)
        ax.set_title(f'{loc.title()} — {ab_label} vs Pre-sorted',
                     fontsize=11, fontweight='bold', color=LOC_C[loc], pad=5)
        ax.spines[['top', 'right']].set_visible(False)
        ax.xaxis.grid(True, alpha=0.2, linestyle='--')
        ax.set_axisbelow(True)

plt.tight_layout()
fig2.savefig(os.path.join(FIG_DIR, 'SupplementaryFigure4.png'),
             dpi=600, bbox_inches='tight', facecolor='white')
plt.close()
print("Saved SupplementaryFigure4.png")
