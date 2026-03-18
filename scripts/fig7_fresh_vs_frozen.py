import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
import warnings; warnings.filterwarnings('ignore')

from shared_utils import (LOC_C, AB_C, LOC_ORDER, AB_LABELS, DATA_DIR, FIG_DIR)

# ── Colour / order constants ───────────────────────────────────────────────────
ORG_ORDER  = [loc.title() for loc in LOC_ORDER]
FRACS_PLOT = ['IgA', 'IgM', 'PreSort']
FRAC_C     = {'IgA': AB_C['A'], 'IgM': AB_C['M'], 'PreSort': AB_C['PRE']}
FRAC_LABEL = {'IgA': AB_LABELS['A'], 'IgM': AB_LABELS['M'], 'PreSort': AB_LABELS['PRE']}

HCMAP = LinearSegmentedColormap.from_list(
    'YlOrDkRd',
    ['#FFF7EC', '#FDD49E', '#FC8D59', '#D7301F', '#7F0000'], N=256)

# ── Load data ─────────────────────────────────────────────────────────────────
bar_df  = pd.read_csv(os.path.join(DATA_DIR, 'lactobacillus_fresh_frozen.csv'))
heat_df = pd.read_csv(os.path.join(DATA_DIR, 'taxa_abundance_fresh_frozen.csv'))


def abbrev(t):
    """Abbreviated italic taxon label for heatmap x-axis."""
    p = t.split()
    if len(p) >= 2:
        return f'$\\it{{{p[0][0]}.~{p[1]}}}$'
    return f'$\\it{{{t}}}$'


# ── Figure layout ─────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(18, 11))
fig.patch.set_facecolor('white')
gs  = fig.add_gridspec(2, 3, hspace=0.25, wspace=0.25, top=0.93, bottom=0.06, left=0.07, right=0.97)
rng = np.random.default_rng(42)
BW  = 0.30

# ══════════════════════════════════════════════════════════════════════════════
# ROW A — Lactobacillus bar charts (lactobacillus_fresh_frozen.csv)
# ══════════════════════════════════════════════════════════════════════════════
for ci, organ in enumerate(ORG_ORDER):
    ax  = fig.add_subplot(gs[0, ci])
    sub = bar_df[bar_df['Organ'] == organ]
    xt, xl = [], []

    for fi, frac in enumerate(FRACS_PLOT):
        fd      = sub[sub['Fraction'] == frac]
        fresh_v = fd[fd['Storage'] == 'Fresh'].groupby('Patient')['Lactobacillus_spp_pct'].sum()
        froz_v  = fd[fd['Storage'] == 'Frozen'].groupby('Patient')['Lactobacillus_spp_pct'].sum()
        col     = FRAC_C[frac]
        n_f, n_z = len(fresh_v), len(froz_v)

        for si, (vals, hatch, alpha) in enumerate([
            (fresh_v, '',    0.85),
            (froz_v,  '///', 0.52),
        ]):
            if not len(vals):
                continue
            xp = fi + (si - 0.5) * BW * 1.08
            ax.bar(xp, vals.mean(), width=BW,
                   color=col, alpha=alpha, hatch=hatch,
                   edgecolor=col if hatch else 'white',
                   linewidth=0.7, zorder=3)
            ax.errorbar(xp, vals.mean(),
                        yerr=vals.std() if len(vals) > 1 else 0,
                        fmt='none', ecolor='#222',
                        elinewidth=1.2, capsize=4, capthick=1.2, zorder=5)
            jit = rng.uniform(-0.06, 0.06, len(vals))
            ax.scatter(xp + jit, vals.values, s=26, color='#444',
                       alpha=0.72, zorder=6, edgecolors='white', lw=0.3)

        xt.append(fi)
        xl.append(f'{FRAC_LABEL[frac]}\n(n={n_f}/{n_z})')

    ax.set_xticks(xt)
    ax.set_xticklabels(xl, fontsize=10, fontweight='bold')
    for tick, frac in zip(ax.get_xticklabels(), FRACS_PLOT):
        tick.set_color(FRAC_C[frac])
    ax.set_ylabel('Lactobacillus spp. abundance (%)', fontsize=9.5)
    ax.set_ylim(0, 125)
    ax.set_yticks([0, 20, 40, 60, 80, 100, 120])
    ax.set_title(organ, fontsize=13, fontweight='bold',
                 color=LOC_C[organ.upper()], pad=5)
    ax.spines[['top', 'right']].set_visible(False)
    ax.yaxis.grid(True, alpha=0.22, linestyle='--', lw=0.5)
    ax.set_axisbelow(True)
    ax.tick_params(axis='y', labelsize=8)

    if ci == 1:
        fp = mpatches.Patch(facecolor='#888', alpha=0.85,
                             edgecolor='white', label='Fresh')
        zp = mpatches.Patch(facecolor='#888', alpha=0.50,
                             hatch='///', edgecolor='#888', label='Frozen')
        ax.legend(handles=[fp, zp], fontsize=9, frameon=True,
                  loc='lower center', ncol=2, framealpha=0.9, edgecolor='#ddd',
                  handlelength=1.8, bbox_to_anchor=(0.5, 0.92))

# ══════════════════════════════════════════════════════════════════════════════
# ROW B — Heatmaps (taxa_abundance_fresh_frozen.csv)
# Columns = top-8 taxa per organ (pre-computed in heatmap file)
# Values  = MeanAbundance_pct
# ══════════════════════════════════════════════════════════════════════════════
VMAX      = 80.0
HEAT_ROWS = [('IgA', 'Fresh'), ('IgA', 'Frozen'), ('IgM', 'Fresh'), ('IgM', 'Frozen'),
             ('PreSort', 'Fresh'), ('PreSort', 'Frozen')]
HEAT_LABS = ['IgA\n(Fresh)', 'IgA\n(Frozen)', 'IgM\n(Fresh)',
             'IgM\n(Frozen)', 'Pre-sort\n(Fresh)', 'Pre-sort\n(Frozen)']

for ci, organ in enumerate(ORG_ORDER):
    ax      = fig.add_subplot(gs[1, ci])
    org_df  = heat_df[heat_df['Organ'] == organ]
    top8    = (org_df.groupby('Taxon')['MeanAbundance_pct']
                     .mean()
                     .nlargest(8)
                     .index.tolist())

    mat = []
    for frac, storage in HEAT_ROWS:
        sg  = org_df[(org_df['Fraction'] == frac) & (org_df['Storage'] == storage)]
        row = []
        for t in top8:
            match = sg[sg['Taxon'] == t]['MeanAbundance_pct']
            row.append(float(match.iloc[0]) if len(match) else 0.0)
        mat.append(row)
    arr = np.array(mat)

    im = ax.imshow(arr, aspect='auto', cmap=HCMAP,
                   vmin=0, vmax=VMAX, interpolation='nearest')

    for sep_y in [1.5, 3.5]:
        ax.axhline(sep_y, color='white', lw=2.2, zorder=3)
    for x in np.arange(-0.5, len(top8), 1):
        ax.axvline(x, color='white', lw=0.7, zorder=2)
    for y in np.arange(-0.5, len(HEAT_ROWS), 1):
        ax.axhline(y, color='white', lw=0.7, zorder=2)

    ax.set_yticks(range(len(HEAT_LABS)))
    ax.set_yticklabels(HEAT_LABS, fontsize=8.5)
    ax.set_xticks(range(len(top8)))
    ax.set_xticklabels([abbrev(t) for t in top8],
                       rotation=42, ha='right', fontsize=8.0)
    ax.tick_params(length=0)
    ax.set_title(organ, fontsize=13, fontweight='bold',
                 color=LOC_C[organ.upper()], pad=5)

    if ci == len(ORG_ORDER) - 1:
        cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cb.set_label('Mean abundance (%)', fontsize=9)
        cb.ax.tick_params(labelsize=8)
        cb.set_ticks([0, 20, 40, 60, 80])

# ── Save ──────────────────────────────────────────────────────────────────────
fig.savefig(os.path.join(FIG_DIR, 'Figure7.png'),
            dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print('Saved Figure7.png')
