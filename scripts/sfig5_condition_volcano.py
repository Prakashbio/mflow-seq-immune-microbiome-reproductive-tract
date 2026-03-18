import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats import mannwhitneyu
import warnings; warnings.filterwarnings('ignore')

from shared_utils import (load_data, LOC_C, LOC_ORDER, FIG_DIR)

otu_rel, meta, tax = load_data()
DATA  = otu_rel.T
taxa  = otu_rel.index.tolist()
PSEUDO = 1e-5

# ─────────────────────────────────────────────────────────────────────────────
# SUPPLEMENTARY: Fold-change volcano Healthy vs Dysbiotic
# ─────────────────────────────────────────────────────────────────────────────
fig2, axes2 = plt.subplots(1, 3, figsize=(20, 8))
fig2.patch.set_facecolor('white')
shared_y_max = 0.0

for oi, loc in enumerate(LOC_ORDER):
    ax    = axes2[oi]
    pre   = meta[(meta['Location'] == loc) & (meta['Antibody'] == 'PRE')]
    h_ids = [s for s in pre[pre['Condition']=='HEALTHY'].index if s in DATA.index]
    d_ids = [s for s in pre[pre['Condition']=='DYSBIOTIC'].index if s in DATA.index]

    rows = []
    for taxon in taxa:
        mu_h  = DATA.loc[h_ids, taxon].mean() if h_ids else 0
        mu_d  = DATA.loc[d_ids, taxon].mean() if d_ids else 0
        lfc   = np.log2((mu_h + PSEUDO) / (mu_d + PSEUDO))

        p_val = 1.0
        if len(h_ids) > 1 and len(d_ids) > 1:
            try:
                _, p_val = mannwhitneyu(
                    DATA.loc[h_ids, taxon].values,
                    DATA.loc[d_ids, taxon].values,
                    alternative='two-sided')
            except: pass
        rows.append({'Taxon': taxon, 'log2FC': lfc, 'p': p_val,
                     'MeanAbund': (mu_h + mu_d) / 2})

    df = pd.DataFrame(rows)
    sig = df['p'] < 0.05
    col_pts = np.where(sig & (df['log2FC'] > 0), '#2196F3',
                np.where(sig & (df['log2FC'] < 0), '#F44336', '#AAAAAA'))

    yvals = -np.log10(df['p'] + 1e-10)
    shared_y_max = max(shared_y_max, float(np.nanmax(yvals)) if len(yvals) else 0.0)

    ax.scatter(df['log2FC'], yvals,
               c=col_pts, s=df['MeanAbund'].clip(0, 1) * 300 + 30,
               alpha=0.80, edgecolors='white', lw=0.5, zorder=5)
    ax.axvline(0, color='#333', lw=1.5, linestyle='--', alpha=0.7)
    ax.axhline(-np.log10(0.05), color='#999', lw=1.2, linestyle=':', alpha=0.8)

    for _, row in df[sig].iterrows():
        ax.annotate(row['Taxon'].split(' ')[1], xy=(row['log2FC'], -np.log10(row['p'] + 1e-10)),
                    xytext=(4, 3), textcoords='offset points', fontsize=7.5, style='italic')

    ax.set_xlabel('log₂(Healthy / Dysbiotic)', fontsize=10)
    ax.set_ylabel('−log₁₀(p)', fontsize=10)
    ax.set_title(loc.title(), fontsize=13, fontweight='bold', color=LOC_C[loc], pad=8)
    ax.spines[['top', 'right']].set_visible(False)
    ax.grid(True, alpha=0.18, linestyle='--')

shared_top = max(shared_y_max * 1.05, -np.log10(0.05) * 1.2)
for ax in axes2:
    ax.set_ylim(0, shared_top)

handles = [mpatches.Patch(color='#2196F3', alpha=0.8, label='Healthy-enriched (p<0.05)'),
           mpatches.Patch(color='#F44336', alpha=0.8, label='Dysbiotic-enriched (p<0.05)'),
           mpatches.Patch(color='#AAAAAA', alpha=0.6, label='Not significant')]
fig2.legend(handles=handles, loc='lower center', ncol=3, fontsize=10,
            frameon=True, bbox_to_anchor=(0.5, -0.03))

plt.tight_layout()
fig2.savefig(os.path.join(FIG_DIR, 'SupplementaryFigure5.png'),
             dpi=600, bbox_inches='tight', facecolor='white')
plt.close()
print("Saved SupplementaryFigure5.png")
