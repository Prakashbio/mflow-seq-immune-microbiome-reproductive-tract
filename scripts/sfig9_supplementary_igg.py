"""
Supplementary Figure 9 – IgG fraction analysis and statistical summary table
=============================================================================
Panel A: IgG vs Pre-sorted comparison (often neglected fraction)
Panel B: IgA + IgM + IgG tri-ternary plot of taxon Ig-coating distribution
Panel C: Pairwise Bray-Curtis distances between all fraction pairs

Statistical summary table is exported as CSV

Output:
  figures/SupplementaryFigure9.png
    figures/SFig9_IgGAnalysis_Table.csv
"""

import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats import kruskal, mannwhitneyu
import warnings; warnings.filterwarnings('ignore')

from shared_utils import (load_data, shannon, hellinger, bray_curtis, pct_label,
                           LOC_C, AB_C, LOC_ORDER, AB_ORDER, AB_LABELS, FIG_DIR)

otu_rel, meta, tax = load_data()
DATA  = otu_rel.T
taxa  = otu_rel.index.tolist()
PSEUDO = 1e-5


# ── Helper: group mean ─────────────────────────────────────────────────────────
def group_mean(loc, ab):
    sids = meta[(meta['Location'] == loc) & (meta['Antibody'] == ab)].index
    sids = [s for s in sids if s in DATA.index]
    if not sids: return pd.Series(0.0, index=taxa)
    return DATA.loc[sids].mean(axis=0)


# ─────────────────────────────────────────────────────────────────────────────
# SFIG 9A: IgG analysis
# ─────────────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(22, 18))
fig.patch.set_facecolor('white')

gs = fig.add_gridspec(3, 3, hspace=0.20, wspace=0.38)

# ── Panels A-C: IgG vs Pre-sorted log2FC per organ ───────────────────────────
for oi, loc in enumerate(LOC_ORDER):
    ax     = fig.add_subplot(gs[0, oi])
    mu_g   = group_mean(loc, 'G')
    mu_pre = group_mean(loc, 'PRE')
    log2fc = np.log2((mu_g.values + PSEUDO) / (mu_pre.values + PSEUDO))
    df     = pd.DataFrame({'Taxon': taxa, 'log2FC': log2fc,
                           'MeanG': mu_g.values * 100}).sort_values('log2FC')
    y      = np.arange(len(df))
    colors_bar = ['#1B9E77' if v > 0 else '#AAAAAA' for v in df['log2FC']]

    ax.barh(y, df['log2FC'], color=colors_bar, alpha=0.75, height=0.65, edgecolor='white')
    ax.scatter(df['log2FC'], y,
               s=df['MeanG'].clip(0, 100) * 3 + 20,
               c=colors_bar, alpha=0.85, zorder=5, edgecolors='white', lw=0.5)
    ax.axvline(0, color='#333', lw=1.5, linestyle='--', alpha=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels(df['Taxon'], fontsize=7.5, style='italic')
    ax.set_xlabel('log₂(IgG / Pre-sorted)', fontsize=9)
    ax.set_title(f'A   {loc.title()} — IgG enrichment', fontsize=11,
                 fontweight='bold', color=LOC_C[loc], loc='left', pad=6)
    ax.spines[['top', 'right']].set_visible(False)
    ax.xaxis.grid(True, alpha=0.2, linestyle='--')
    ax.set_axisbelow(True)

# ── Panel B: Ternary-style IgA/IgM/IgG scatter per taxon ─────────────────────
# For each taxon: compute proportional IgA / (IgA+IgM+IgG) as x-axis
# and IgM / (IgA+IgM+IgG) as y-axis; IgG = remainder
ax_tern = fig.add_subplot(gs[1, :])

all_loc_data = []
for loc in LOC_ORDER:
    for ti, taxon in enumerate(taxa):
        mu_a = group_mean(loc, 'A').get(taxon, 0) * 100
        mu_m = group_mean(loc, 'M').get(taxon, 0) * 100
        mu_g = group_mean(loc, 'G').get(taxon, 0) * 100
        tot  = mu_a + mu_m + mu_g
        if tot < 0.3: continue
        all_loc_data.append({
            'Taxon': taxon, 'Location': loc,
            'IgA_prop': mu_a / tot, 'IgM_prop': mu_m / tot, 'IgG_prop': mu_g / tot,
            'TotalAbund': tot
        })

df_tern = pd.DataFrame(all_loc_data)
if not df_tern.empty:
    sc = ax_tern.scatter(df_tern['IgA_prop'], df_tern['IgM_prop'],
                         s=df_tern['TotalAbund'].clip(0, 50) * 8 + 30,
                         c=[LOC_C[l] for l in df_tern['Location']],
                         alpha=0.72, edgecolors='white', lw=0.8, zorder=5)
    label_offsets = {
        'VAGINA': (10, 8),
        'CERVIX': (10, -10),
        'ENDOMETRIUM': (-12, 8),
    }

    # Diagonal reference lines
    ax_tern.axhline(0.33, color='#ccc', lw=1, linestyle=':')
    ax_tern.axvline(0.33, color='#ccc', lw=1, linestyle=':')
    ax_tern.plot([0, 1], [1, 0], 'k--', lw=0.8, alpha=0.25)

    # Label key taxa
    for taxon in ['Lactobacillus iners', 'Gardnerella vaginalis', 'Sneathia vaginalis',
                  'Lactobacillus crispatus', 'Parvimonas micra']:
        sub_t = df_tern[df_tern['Taxon'] == taxon]
        if sub_t.empty: continue
        for _, row in sub_t.iterrows():
            dx, dy = label_offsets.get(row['Location'], (8, 6))
            ax_tern.annotate(taxon.split(' ')[1], xy=(row['IgA_prop'], row['IgM_prop']),
                             xytext=(dx, dy), textcoords='offset points',
                             fontsize=8, style='italic', zorder=9,
                             bbox=dict(boxstyle='round,pad=0.15', fc='white', ec='none', alpha=0.80),
                             arrowprops=dict(arrowstyle='-', color='#999', lw=0.6))

    ax_tern.set_xlabel('IgA proportion (IgA / IgA+IgM+IgG)', fontsize=11)
    ax_tern.set_ylabel('IgM proportion (IgM / IgA+IgM+IgG)', fontsize=11)
    ax_tern.set_title(
        'B   Proportional Ig-coating space — each point = one taxon per organ (top-right = IgA-dominated; top-left = IgM-dominated; bottom = IgG-dominated)',
        fontsize=11, fontweight='bold', loc='left', pad=6)
    ax_tern.spines[['top', 'right']].set_visible(False)
    ax_tern.grid(True, alpha=0.18, linestyle='--')
    ax_tern.set_xlim(-0.05, 1.05)
    ax_tern.set_ylim(-0.05, 1.05)

    # Annotations for regions
    ax_tern.text(0.85, 0.1, 'IgA\ndominant', color='#E66101',
                 fontsize=10, fontweight='bold', ha='center')
    ax_tern.text(0.1, 0.85, 'IgM\ndominant', color='#5E3C99',
                 fontsize=10, fontweight='bold', ha='center')
    ax_tern.text(0.1, 0.1, 'IgG\ndominant', color='#1B9E77',
                 fontsize=10, fontweight='bold', ha='center')

    loc_handles = [mpatches.Patch(color=LOC_C[l], label=l.title()) for l in LOC_ORDER]
    ax_tern.legend(handles=loc_handles, fontsize=10, frameon=True, loc='upper right')

# ── Export pairwise comparison table as CSV ──────────────────────────────────
div_rows = [{'SampleID': s, 'Shannon': shannon(DATA.loc[s].values)}
            for s in DATA.index]
div = pd.DataFrame(div_rows).set_index('SampleID').join(meta)

pair_list = [('A', 'M'), ('A', 'G'), ('A', 'PRE'), ('M', 'G'), ('M', 'PRE'), ('G', 'PRE')]
table_data = []

for loc in LOC_ORDER:
    for (ab1, ab2) in pair_list:
        g1 = div[(div['Location'] == loc) & (div['Antibody'] == ab1)]['Shannon'].dropna()
        g2 = div[(div['Location'] == loc) & (div['Antibody'] == ab2)]['Shannon'].dropna()
        if len(g1) > 1 and len(g2) > 1:
            try:
                U, p = mannwhitneyu(g1, g2, alternative='two-sided')
                table_data.append([loc.title(),
                                   f"{AB_LABELS[ab1]} vs {AB_LABELS[ab2]}",
                                   f"n={len(g1)}", f"n={len(g2)}",
                                   f"{U:.0f}", f"{p:.4f}", pct_label(p)])
            except: pass

if table_data:
    col_labels_tab = ['Site', 'Comparison', 'n (grp1)', 'n (grp2)', 'U stat', 'p-value', 'Sig']
    table_df = pd.DataFrame(table_data, columns=col_labels_tab)
    table_df.to_csv(os.path.join(FIG_DIR, 'IgGAnalysis_Table.csv'), index=False)


# ── Panel C: Pairwise Bray-Curtis distances between fraction pairs ───────────
PAIR_COMPS = [('A', 'M'), ('A', 'G'), ('A', 'PRE'), ('M', 'G'), ('M', 'PRE'), ('G', 'PRE')]
PAIR_COLS  = ['#D6604D', '#74C476', '#4393C3', '#9E9AC8', '#FDAE6B', '#AAAAAA']

for oi, loc in enumerate(LOC_ORDER):
    ax = fig.add_subplot(gs[2, oi])
    box_data, box_labels, box_colors = [], [], []

    for (ab1, ab2), col in zip(PAIR_COMPS, PAIR_COLS):
        sids1 = [s for s in meta[(meta['Location'] == loc) & (meta['Antibody'] == ab1)].index
                 if s in DATA.index]
        sids2 = [s for s in meta[(meta['Location'] == loc) & (meta['Antibody'] == ab2)].index
                 if s in DATA.index]
        if not sids1 or not sids2:
            continue

        X1 = hellinger(otu_rel[sids1].T.values)
        X2 = hellinger(otu_rel[sids2].T.values)
        from scipy.spatial.distance import cdist
        D_cross = cdist(X1, X2, metric='braycurtis').flatten()
        box_data.append(D_cross)
        box_labels.append(f"{AB_LABELS[ab1]}\nvs\n{AB_LABELS[ab2]}")
        box_colors.append(col)

    if not box_data:
        continue

    bp = ax.boxplot(box_data, patch_artist=True, showfliers=True,
                    medianprops={'color': 'white', 'lw': 2.5},
                    flierprops={'marker': 'o', 'markersize': 3, 'alpha': 0.5})
    for patch, col in zip(bp['boxes'], box_colors):
        patch.set_facecolor(col)
        patch.set_alpha(0.6)
    for i, (whisker, cap) in enumerate(zip(bp['whiskers'], bp['caps'])):
        c = box_colors[i // 2]
        whisker.set_color(c)
        cap.set_color(c)

    ax.set_xticks(range(1, len(box_labels) + 1))
    ax.set_xticklabels(box_labels, fontsize=8)
    ax.set_ylabel('Bray-Curtis dissimilarity', fontsize=10)
    ax.set_title(f"{'C   ' if oi == 0 else ''}{loc.title()} — pairwise fraction distances",
                 fontsize=11 if oi == 0 else 13,
                 fontweight='bold', color=LOC_C[loc], loc='left' if oi == 0 else 'center', pad=8)
    ax.spines[['top', 'right']].set_visible(False)
    ax.yaxis.grid(True, alpha=0.22, linestyle='--')
    ax.set_axisbelow(True)
    ax.set_ylim(0, 1)

plt.tight_layout()
fig.savefig(os.path.join(FIG_DIR, 'SupplementaryFigure9.png'),
            dpi=600, bbox_inches='tight', facecolor='white')
plt.close()
print("Saved SupplementaryFigure9.png")
