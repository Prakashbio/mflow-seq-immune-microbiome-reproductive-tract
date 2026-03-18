import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import warnings; warnings.filterwarnings('ignore')

from shared_utils import load_data, LOC_C, AB_C, LOC_ORDER, FIG_DIR


otu_rel, meta, tax = load_data()
DATA   = otu_rel.T
TAXA   = otu_rel.index.tolist()
PSEUDO = 1e-5


def group_mean(loc, ab):
    sids = meta[(meta.Location == loc) & (meta.Antibody == ab)].index
    sids = [s for s in sids if s in DATA.index]
    if not sids:
        return pd.Series(0.0, index=TAXA), 0
    return DATA.loc[sids].mean(axis=0) * 100, len(sids)


def abbrev(name):
    name_clean = name.replace('†', '').strip()
    parts = name_clean.split()
    if len(parts) >= 2:
        return f"$\\it{{{parts[0][0]}.~{'~'.join(parts[1:])}}}$"
    return f"$\\it{{{name_clean}}}$"



fig = plt.figure(figsize=(26, 12))
fig.patch.set_facecolor('white')

gs = fig.add_gridspec(1, 3, hspace=0.15, wspace=0.25)

# ── Draw each panel
for ci, loc in enumerate(LOC_ORDER):
    ax = fig.add_subplot(gs[0, ci])

    mu_a, n_a = group_mean(loc, 'A')
    mu_m, n_m = group_mean(loc, 'M')
    mu_p, n_p = group_mean(loc, 'PRE')

    # Differential abundance: IgA - IgM
    diff = mu_a - mu_m

    # Build dataframe
    df = pd.DataFrame({
        'Taxon': TAXA,
        'IgA': mu_a,
        'IgM': mu_m,
        'Pre': mu_p,
        'Diff': diff,
        'Cont': False
    })

    # Add contaminant for Cervix (transparency)
    if loc == 'CERVIX':
        df = pd.concat([df, pd.DataFrame([{
            'Taxon': 'Mageibacillus indolicus',
            'IgA': 0.3,
            'IgM': 60.3,
            'Pre': 6.0,
            'Diff': -60.0,
            'Cont': True
        }])], ignore_index=True)

    df = df.sort_values('Diff').reset_index(drop=True)
    n = len(df)
    y = np.arange(n)

    # Colors
    colors = [('#AAAAAA' if r.Cont else
               (AB_C['A'] if r.Diff >= 0 else AB_C['M']))
              for _, r in df.iterrows()]

    # ── Horizontal bars ───────────────────────────────────────────────────────
    ax.barh(y, df['Diff'], color=colors, alpha=0.40, height=0.50,
            edgecolor='white', linewidth=0.3, zorder=2)

    # ── Dots at bar tips ──────────────────────────────────────────────────────
    # Size = pre-sorted abundance (clipped at 80%)
    dot_size = np.clip(df['Pre'].values, 0, 80) * 6 + 20

    norm_mask = ~df['Cont'].values
    if norm_mask.any():
        ax.scatter(df['Diff'].values[norm_mask], y[norm_mask],
                   c=[c for c, m in zip(colors, norm_mask) if m],
                   s=dot_size[norm_mask], zorder=5,
                   edgecolors='white', linewidths=0.6, alpha=0.90)

    cont_mask = df['Cont'].values
    if cont_mask.any():
        ax.scatter(df['Diff'].values[cont_mask], y[cont_mask],
                   c='#AAAAAA', s=dot_size[cont_mask], zorder=5,
                   edgecolors='#666', linewidths=1.5, alpha=0.75)

    # ── Zero reference line ───────────────────────────────────────────────────
    ax.axvline(0, color='#555', lw=1.2, linestyle='--', alpha=0.55, zorder=3)

    # ── Grid ──────────────────────────────────────────────────────────────────
    ax.xaxis.grid(True, alpha=0.18, linestyle='--', linewidth=0.6, zorder=0)
    ax.set_axisbelow(True)

    # ── Y-axis: taxon names ───────────────────────────────────────────────────
    ax.set_yticks(y)
    ylabels = []
    for _, r in df.iterrows():
        lbl = abbrev(r['Taxon'])
        if r['Cont']:
            lbl = lbl + '†'
        ylabels.append(lbl)
    ax.set_yticklabels(ylabels, fontsize=9, style='normal')
    for tick, cont in zip(ax.get_yticklabels(), df['Cont']):
        if cont:
            tick.set_color('#AAAAAA')
    ax.tick_params(axis='y', length=0, pad=2)

    # ── X-axis ────────────────────────────────────────────────────────────────
    ax.set_xlabel('IgA abundance (%) − IgM abundance (%)', fontsize=8.5, labelpad=4)

    # ── Corner enrichment labels with sample counts ───────────────────────────
    xl = ax.get_xlim()
    ax.text(xl[1] * 0.92, -1.2, f'IgA-enriched → (n={n_a})',
            color=AB_C['A'], fontsize=8, fontweight='bold',
            ha='right', va='bottom')
    ax.text(0.02, 0.98, f'← IgM-enriched (n={n_m})',
            transform=ax.transAxes,
            color=AB_C['M'], fontsize=8, fontweight='bold',
            ha='left', va='top')

    # ── Title ─────────────────────────────────────────────────────────────────
    ax.set_title(loc.title(), fontsize=11, fontweight='bold',
                 color=LOC_C[loc], pad=5)

    # ── Spines ────────────────────────────────────────────────────────────────
    ax.spines[['top', 'right']].set_visible(False)
    ax.spines[['left', 'bottom']].set_linewidth(0.8)
    ax.set_ylim(-1.5, n + 0.3)

# ── Figure-level legend ───────────────────────────────────────────────────────
legend_handles = [
    mpatches.Patch(facecolor=AB_C['A'], alpha=0.75,
                   edgecolor='none', label='IgA-enriched'),
    mpatches.Patch(facecolor=AB_C['M'], alpha=0.75,
                   edgecolor='none', label='IgM-enriched'),
    mpatches.Patch(facecolor='#AAAAAA', alpha=0.65,
                   edgecolor='#666', linewidth=1.5,
                   label='Excluded contaminant (†)'),
]
# Dot-size key
for pct, lbl in [(2, '2%'), (20, '20%'), (50, '50%')]:
    sz = np.clip(pct, 0, 80) * 6 + 20
    legend_handles.append(
        Line2D([0], [0], marker='o', color='none',
               markerfacecolor='#888888', markeredgecolor='white',
               markeredgewidth=0.6, markersize=np.sqrt(sz),
               label=f'Pre-sorted: {lbl}')
    )

fig.legend(handles=legend_handles,
           loc='lower center', ncol=6,
           fontsize=8, frameon=True, framealpha=0.95, edgecolor='#ddd',
           bbox_to_anchor=(0.5, 0.05),
           handletextpad=0.5, columnspacing=1.2)

plt.tight_layout()

# ── Save ──────────────────────────────────────────────────────────────────────
fig.savefig(os.path.join(FIG_DIR, 'Figure8.png'),
            dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print('Saved Figure8.png')
