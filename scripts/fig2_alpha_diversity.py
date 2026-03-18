import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from scipy.stats import mannwhitneyu
import warnings; warnings.filterwarnings('ignore')

from shared_utils import (load_data, shannon, simpson, richness, evenness,
                           LOC_C, AB_C, COND_C, LOC_ORDER, AB_ORDER, AB_LABELS, pct_label, FIG_DIR)

otu_rel, meta, tax = load_data()
DATA = otu_rel.T  # samples × taxa

# ── Compute per-sample diversity ──────────────────────────────────────────────
div_rows = []
for sid, row in DATA.iterrows():
    counts = row.values
    div_rows.append({
        'SampleID': sid,
        'Shannon':  shannon(counts),
        'Simpson':  simpson(counts),
        'Richness': richness(counts),
        'Evenness': evenness(counts),
    })
div = pd.DataFrame(div_rows).set_index('SampleID')
div = div.join(meta)

METRICS      = ['Shannon', 'Simpson', 'Richness', 'Evenness']
METRIC_YLABS = ["Shannon H'", "Simpson (1-D)", "Taxon richness", "Pielou's evenness J'"]


def draw_alpha_panels(div_sub, strat_col, strat_vals, strat_colors, fname):
    """Generic function to draw alpha diversity figure."""
    fig, axes = plt.subplots(len(METRICS), 3, figsize=(20, 16))
    fig.patch.set_facecolor('white')

    for mi, (metric, ylab) in enumerate(zip(METRICS, METRIC_YLABS)):
        for oi, loc in enumerate(LOC_ORDER):
            ax      = axes[mi, oi]
            sub     = div_sub[div_sub['Location'] == loc]
            n_strat = len(strat_vals)
            n_ab    = len(AB_ORDER)
            bw      = 0.35

            all_vals = []
            for ai, ab in enumerate(AB_ORDER):
                x_base = ai * (n_strat + 0.5)
                for si, sv in enumerate(strat_vals):
                    grp = sub[(sub['Antibody'] == ab) & (sub[strat_col] == sv)][metric].dropna()
                    if len(grp) == 0:
                        continue
                    xpos = x_base + si * bw
                    bp = ax.boxplot(grp, positions=[xpos], widths=bw * 0.75,
                                    patch_artist=True, showfliers=False,
                                    medianprops={'color': 'white', 'lw': 2, 'zorder': 10},
                                    whiskerprops={'color': strat_colors[sv], 'lw': 1.5},
                                    capprops=  {'color': strat_colors[sv], 'lw': 1.5},
                                    boxprops=  {'facecolor': strat_colors[sv], 'alpha': 0.55,
                                                 'edgecolor': strat_colors[sv], 'lw': 1.2})
                    np.random.seed(ai * 10 + si)
                    jit = np.random.normal(0, 0.04, len(grp))
                    ax.scatter(xpos + jit, grp, c=strat_colors[sv], s=35,
                               alpha=0.8, zorder=8, edgecolors='white', lw=0.3)
                    all_vals.extend(grp.tolist())

            # Pairwise significance per fraction between the two strata
            if n_strat == 2 and all_vals:
                y_span = max(max(all_vals) - min(all_vals), 1e-6)
                y_pad = 0.06 * y_span
                y_lift = 0.02 * y_span
                y_text = 0.015 * y_span
                y_top = max(all_vals)

                for ai, ab in enumerate(AB_ORDER):
                    x_base = ai * (n_strat + 0.5)
                    x1 = x_base
                    x2 = x_base + bw

                    g1 = sub[(sub['Antibody'] == ab) &
                             (sub[strat_col] == strat_vals[0])][metric].dropna()
                    g2 = sub[(sub['Antibody'] == ab) &
                             (sub[strat_col] == strat_vals[1])][metric].dropna()

                    if len(g1) > 0 and len(g2) > 0:
                        try:
                            _, p = mannwhitneyu(g1, g2, alternative='two-sided')
                            sig_lbl = pct_label(p)
                        except Exception:
                            sig_lbl = 'NA'
                        group_max = max(g1.max(), g2.max())
                    else:
                        sig_lbl = 'NA'
                        group_max = y_top

                    y = group_max + y_pad
                    ax.plot([x1, x1, x2, x2],
                            [y, y + y_lift, y + y_lift, y],
                            color='black', lw=0.9, zorder=11)
                    ax.text((x1 + x2) / 2, y + y_lift + y_text, sig_lbl,
                            ha='center', va='bottom', fontsize=8,
                            fontweight='bold', color='black')
                    y_top = max(y_top, y + y_lift + 2 * y_text)

                ymin, ymax = ax.get_ylim()
                if y_top > ymax:
                    ax.set_ylim(ymin, y_top + y_pad)

            # x-ticks
            tick_pos = [ai * (n_strat + 0.5) + (n_strat - 1) * bw / 2
                        for ai in range(n_ab)]
            ax.set_xticks(tick_pos)
            ax.set_xticklabels([AB_LABELS[ab] for ab in AB_ORDER],
                               fontsize=10, fontweight='bold')
            for tick, ab in zip(ax.get_xticklabels(), AB_ORDER):
                tick.set_color(AB_C[ab])

            ax.set_ylabel(ylab, fontsize=11)
            if mi == 0:
                ax.set_title(loc.title(), fontsize=13, fontweight='bold',
                             color=LOC_C[loc], pad=8)
            ax.spines[['top', 'right']].set_visible(False)
            ax.yaxis.grid(True, alpha=0.22, linestyle='--')
            ax.set_axisbelow(True)

    # Legend
    handles = [mpatches.Patch(color=strat_colors[sv], alpha=0.75,
                              label=sv.title()) for sv in strat_vals]
    sig_handles = [
        Line2D([], [], color='none', label='Significance:'),
        Line2D([], [], color='none', label='ns (p >= 0.05)'),
        Line2D([], [], color='none', label='* (p < 0.05)'),
    ]
    handles.extend(sig_handles)
    fig.legend(handles=handles, loc='lower center', ncol=len(handles),
               fontsize=11, frameon=True, bbox_to_anchor=(0.5, -0.02))

    plt.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, fname), dpi=600,
                bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Saved {fname}")


# ── Main figure: stratify by Condition ───────────────────────────────────────
draw_alpha_panels(
    div, 'Condition', ['HEALTHY', 'DYSBIOTIC'], COND_C,
    'Figure2.png',
)
