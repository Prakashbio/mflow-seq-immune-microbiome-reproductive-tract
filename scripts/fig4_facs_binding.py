import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from collections import defaultdict
from scipy.stats import mannwhitneyu
import warnings; warnings.filterwarnings('ignore')

from shared_utils import (load_data, pct_label, LOC_C, AB_C,
                          LOC_ORDER, AB_LABELS, DATA_DIR, FIG_DIR)

# ── Colour map: Ig fraction label → colour ────────────────────────────────────
FRAC_ORDER  = ['IgA', 'IgM', 'IgG']
FRAC_COLORS = {'IgA': AB_C['A'], 'IgM': AB_C['M'], 'IgG': AB_C['G']}

# ── Load FACS cell counts ──────────────────────────────────────────────────────
facs_df = pd.read_csv(os.path.join(DATA_DIR, 'facs_binding_events.csv'))
cells   = defaultdict(list)
for _, row in facs_df.iterrows():
    try:
        val = float(row['FACSEvents'])
        if val > 0:
            cells[(row['Location'], row['Fraction'])].append(val)
    except (ValueError, TypeError):
        pass

# ── Load sequencing reads ──────────────────────────────────────────────────────
otu_raw      = pd.read_csv(os.path.join(DATA_DIR, 'otu_table.csv'), index_col=0)
_, meta, _   = load_data()

reads = defaultdict(list)
for sid in meta.index:
    ab = meta.loc[sid, 'Antibody']
    if ab not in ('A', 'M', 'G'):
        continue
    frac = AB_LABELS[ab]
    loc  = meta.loc[sid, 'Location']
    val  = float(otu_raw[sid].sum()) if sid in otu_raw.columns else 0.0
    if val > 0:
        reads[(loc, frac)].append(val)

# ── Helpers ───────────────────────────────────────────────────────────────────
def wilcox_p(a, b):
    a = [x for x in a if x > 0]
    b = [x for x in b if x > 0]
    if len(a) < 3 or len(b) < 3:
        return 1.0
    _, p = mannwhitneyu(a, b, alternative='two-sided')
    return p


def violin_panel(ax, data_dict, site, ylabel, title_str, add_pval=True):
    """Draw a violin + jitter + median panel for one anatomical site."""
    positions = [1, 2, 3]
    vals_list = [data_dict.get((site, frac), []) for frac in FRAC_ORDER]
    log_vals  = [
        np.log10(np.array([x for x in v if x > 0])) if any(x > 0 for x in v)
        else np.array([0.0])
        for v in vals_list
    ]

    # ── Violin bodies ─────────────────────────────────────────────────────────
    parts = ax.violinplot(log_vals, positions=positions,
                          showmedians=False, showextrema=False, widths=0.75)
    for i, pc in enumerate(parts['bodies']):
        pc.set_facecolor(FRAC_COLORS[FRAC_ORDER[i]])
        pc.set_alpha(0.50)
        pc.set_edgecolor('none')
        pc.set_zorder(2)

    # ── Median bar ────────────────────────────────────────────────────────────
    for i, lv in enumerate(log_vals):
        if len(lv) == 0:
            continue
        med = np.median(lv)
        ax.hlines(med, positions[i] - 0.22, positions[i] + 0.22,
                  colors='black', linewidths=2.8, zorder=5)

    # ── Jittered dots ─────────────────────────────────────────────────────────
    rng = np.random.default_rng(42)
    for i, (frac, lv) in enumerate(zip(FRAC_ORDER, log_vals)):
        if len(lv) == 0:
            continue
        jitter = rng.uniform(-0.10, 0.10, len(lv))
        ax.scatter(positions[i] + jitter, lv,
                   color=FRAC_COLORS[frac], s=28, alpha=0.90,
                   zorder=6, edgecolors='white', linewidths=0.5)

    # ── IgA vs IgM significance bracket ──────────────────────────────────────
    if add_pval and len(vals_list[0]) >= 3 and len(vals_list[1]) >= 3:
        p      = wilcox_p(vals_list[0], vals_list[1])
        lbl    = f'IgA vs IgM: {pct_label(p)}'
        all_lv = np.concatenate([lv for lv in log_vals if len(lv)])
        y_br   = all_lv.max() + 0.30
        ax.plot([1, 1, 2, 2], [y_br - 0.08, y_br, y_br, y_br - 0.08],
                lw=0.9, color='black', zorder=7)
        ax.text(1.5, y_br + 0.04, lbl, ha='center', va='bottom',
                fontsize=7.5, color='black')

    # ── Y-axis log ticks ──────────────────────────────────────────────────────
    all_lv = np.concatenate([lv for lv in log_vals if len(lv)])
    ymin = np.floor(all_lv.min()) - 0.2
    ymax = np.ceil(all_lv.max()) + 0.6
    ax.set_ylim(ymin, ymax)
    yticks = np.arange(np.ceil(ymin), np.floor(ymax) + 1, 1.0)
    ax.set_yticks(yticks)
    ax.set_yticklabels([f'$10^{{{int(t)}}}$' for t in yticks], fontsize=8)
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.tick_params(axis='y', which='minor', length=2)

    # ── X-axis labels with per-fraction n counts ──────────────────────────────
    xlabels = [f'{frac}\n(n={len(v)})' for frac, v in zip(FRAC_ORDER, vals_list)]
    ax.set_xticks(positions)
    ax.set_xticklabels(xlabels, fontsize=10, fontweight='bold')
    for tick, frac in zip(ax.get_xticklabels(), FRAC_ORDER):
        tick.set_color(FRAC_COLORS[frac])

    ax.set_ylabel(ylabel, fontsize=9)
    ax.set_title(title_str, fontsize=12, fontweight='bold',
                 color=LOC_C.get(site, '#333333'), pad=6)
    ax.spines[['top', 'right']].set_visible(False)
    ax.spines[['left', 'bottom']].set_linewidth(0.8)
    ax.tick_params(axis='both', which='major', labelsize=8)
    ax.yaxis.grid(True, alpha=0.20, linestyle='--', linewidth=0.5)
    ax.set_axisbelow(True)


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE: FACS cell counts (row 0) and sequencing reads (row 1) per site
# ─────────────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 3, figsize=(13, 9),
                         gridspec_kw={'hspace': 0.15, 'wspace': 0.25})
fig.patch.set_facecolor('white')

for col, site in enumerate(LOC_ORDER):
    violin_panel(axes[0, col], cells, site,
                 ylabel='FACS-sorted bacterial events',
                 title_str=site.title(), add_pval=True)
    violin_panel(axes[1, col], reads, site,
                 ylabel='Microbial sequencing reads',
                 title_str='', add_pval=True)

plt.subplots_adjust(top=0.93, bottom=0.07)
fig.savefig(os.path.join(FIG_DIR, 'Figure4.png'),
            dpi=600, bbox_inches='tight', facecolor='white')
plt.close()
print('Saved Figure4.png')
