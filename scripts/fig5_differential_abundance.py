import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import warnings; warnings.filterwarnings('ignore')
from shared_utils import (load_data, LOC_C, LOC_ORDER, FIG_DIR)

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
# MAIN FIGURE
# ─────────────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(22, 20))
fig.patch.set_facecolor('white')

gs = fig.add_gridspec(3, 3, hspace=0.15, wspace=0.35)

# ── Panel A: Heatmap ─────────────────────────────────────────────────────────
ax_heat = fig.add_subplot(gs[0, :])

# Build matrix: taxa × (loc×fraction) for A and M
col_labels = []
heat_data  = []
ab_label_lower = {'A': 'IgA', 'M': 'IgM', 'PRE': 'Pre-sorted'}
for loc in LOC_ORDER:
    for ab in ['A', 'M', 'PRE']:
        col_labels.append(f"{ab_label_lower[ab]}\n{loc.title()}")
        heat_data.append(group_mean(loc, ab).values * 100)

heat_mat = np.array(heat_data).T   # taxa × groups

# Sort taxa by total abundance descending
sort_idx = heat_mat.sum(axis=1).argsort()[::-1]
heat_mat_s = heat_mat[sort_idx]
taxa_s     = [taxa[i] for i in sort_idx]

im = ax_heat.imshow(heat_mat_s, aspect='auto', cmap='YlOrRd',
                    interpolation='nearest', vmin=0, vmax=heat_mat_s.max() * 0.7)
ax_heat.set_xticks(range(len(col_labels)))
ax_heat.set_xticklabels(col_labels, fontsize=9)
ax_heat.set_yticks(range(len(taxa_s)))
ax_heat.set_yticklabels(taxa_s, fontsize=8, style='italic')
ax_heat.set_title('A   Mean relative abundance (%) by anatomical site and Ig fraction',
                  fontsize=11, fontweight='bold', loc='left', pad=6)

# Organ dividers
for xi in [2.5, 5.5]:
    ax_heat.axvline(xi, color='white', lw=3)

plt.colorbar(im, ax=ax_heat, label='Mean rel. abundance (%)',
             fraction=0.015, pad=0.01)


# ── Panels B: Log2FC lollipop per location ────────────────────────────────────
for oi, loc in enumerate(LOC_ORDER):
    ax = fig.add_subplot(gs[1, oi])

    mu_a   = group_mean(loc, 'A')
    mu_m   = group_mean(loc, 'M')
    log2fc = np.log2((mu_a.values + PSEUDO) / (mu_m.values + PSEUDO))

    df_lfc = pd.DataFrame({'Taxon': taxa, 'log2FC': log2fc,
                            'MeanAbund': (mu_a.values + mu_m.values) / 2})
    df_lfc = df_lfc.sort_values('log2FC')

    colors = ['#E66101' if v > 0 else '#5E3C99' for v in df_lfc['log2FC']]
    y = np.arange(len(df_lfc))
    ax.barh(y, df_lfc['log2FC'], color=colors, alpha=0.70, height=0.65,
            edgecolor='white')
    ax.scatter(df_lfc['log2FC'], y,
               c=colors,
               s=df_lfc['MeanAbund'].clip(0, 1) * 400 + 20,
               zorder=5, edgecolors='white', lw=0.6)
    ax.axvline(0, color='#333', lw=1.5, linestyle='--', alpha=0.7)

    ax.set_yticks(y)
    ax.set_yticklabels(df_lfc['Taxon'], fontsize=7.5, style='italic')
    ax.set_xlabel('log₂(IgA / IgM)', fontsize=9)
    ax.set_title(f'B   {loc.title()}', fontsize=11, fontweight='bold',
                 color=LOC_C[loc], loc='left', pad=6)
    ax.spines[['top', 'right']].set_visible(False)
    ax.xaxis.grid(True, alpha=0.2, linestyle='--')
    ax.set_axisbelow(True)

    xl = ax.get_xlim()
    ax.text(xl[1] * 0.92, -1.2, 'IgA-\nenriched', color='#E66101',
            fontsize=8, fontweight='bold', ha='right', va='bottom')
    ax.text(0.02, 0.98, 'IgM-\nenriched', transform=ax.transAxes,
            color='#5E3C99', fontsize=8, fontweight='bold',
            ha='left', va='top')


# ── Panel C: Bubble chart IgA vs IgM (all locations combined) ────────────────
ax_bubble = fig.add_subplot(gs[2, :])

# Per-taxon per-location
for oi, loc in enumerate(LOC_ORDER):
    mu_a = group_mean(loc, 'A').values * 100
    mu_m = group_mean(loc, 'M').values * 100
    mu_p = group_mean(loc, 'PRE').values * 100

    def bubble_color(a, m):
        if a > m * 1.3:  return '#E66101'
        if m > a * 1.3:  return '#5E3C99'
        return '#888888'

    for ti, t in enumerate(taxa):
        if mu_a[ti] + mu_m[ti] < 0.5:
            continue
        col  = bubble_color(mu_a[ti], mu_m[ti])
        size = max(mu_p[ti] * 15 + 20, 20)
        ax_bubble.scatter(mu_a[ti], mu_m[ti], s=size, c=col, alpha=0.65,
                          edgecolors=LOC_C[loc], linewidths=1.2, zorder=5)
        if mu_a[ti] + mu_m[ti] > 5:
            ax_bubble.annotate(t, xy=(mu_a[ti], mu_m[ti]),
                               xytext=(5, 3), textcoords='offset points',
                               fontsize=7, style='italic', color='#222',
                               arrowprops=dict(arrowstyle='-', color='#ccc', lw=0.5))

maxv = max(ax_bubble.get_xlim()[1], ax_bubble.get_ylim()[1]) * 1.05
ax_bubble.plot([0, maxv], [0, maxv], 'k--', lw=1.2, alpha=0.35)
ax_bubble.set_xlim(-0.5, maxv)
ax_bubble.set_ylim(-0.5, maxv)
ax_bubble.set_xlabel('Mean relative abundance in IgA fraction (%)', fontsize=11)
ax_bubble.set_ylabel('Mean relative abundance in IgM fraction (%)', fontsize=11)
ax_bubble.set_title('C   IgA vs IgM taxon enrichment bubble chart  '
                    '(bubble size = pre-sorted abundance; border = anatomical site)',
                    fontsize=11, fontweight='bold', loc='left', pad=6)
ax_bubble.spines[['top', 'right']].set_visible(False)
ax_bubble.grid(True, alpha=0.18, linestyle='--')

# Legends
loc_handles = [mpatches.Patch(facecolor='white', edgecolor=LOC_C[l],
                              label=l.title(), lw=2) for l in LOC_ORDER]
dir_handles = [mpatches.Patch(color='#E66101', alpha=0.75, label='IgA-enriched'),
               mpatches.Patch(color='#5E3C99', alpha=0.75, label='IgM-enriched'),
               mpatches.Patch(color='#888888', alpha=0.75, label='Mixed')]
ax_bubble.legend(handles=loc_handles + dir_handles, fontsize=9, frameon=True,
                 loc='upper left', ncol=3)

plt.tight_layout()
fig.savefig(os.path.join(FIG_DIR, 'Figure5.png'),
            dpi=600, bbox_inches='tight', facecolor='white')
plt.close()
print("Saved Figure5.png")
