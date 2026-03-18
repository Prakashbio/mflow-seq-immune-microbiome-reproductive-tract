import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import warnings; warnings.filterwarnings('ignore')

from shared_utils import load_data, LOC_C, AB_C, LOC_ORDER, AB_ORDER, AB_LABELS, FIG_DIR

# ── Load ──────────────────────────────────────────────────────────────────────
otu_rel, meta, tax = load_data()
DATA = otu_rel.T  # (samples × taxa)

# ── Colour palette for 41 taxa ────────────────────────────────────────────────
TAXA_COLORS = [
    '#1F78B4','#A6CEE3','#33A02C','#B2DF8A','#E31A1C','#FB9A99',
    '#FF7F00','#FDBF6F','#6A3D9A','#CAB2D6','#B15928','#FFFF99',
    '#D95F02','#7570B3','#E7298A','#66A61E','#E6AB02','#A6761D',
    '#2166AC','#4DAC26','#D7191C','#636363','#1B9E77','#D6604D',
    '#4393C3','#F4A582','#92C5DE','#F7F7F7','#BABABA','#404040',
    '#8DD3C7','#FFFFB3','#BEBADA','#FB8072','#80B1D3','#FDB462',
    '#B3DE69','#FCCDE5','#D9D9D9','#BC80BD','#CCEBC5',
]

taxa_list = otu_rel.index.tolist()
taxa_pal  = {t: TAXA_COLORS[i % len(TAXA_COLORS)] for i, t in enumerate(taxa_list)}

# ─────────────────────────────────────────────────────────────────────────────
# MAIN FIGURE: Species composition
# ─────────────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(22, 9))
fig.patch.set_facecolor('white')

for oi, loc in enumerate(LOC_ORDER):
    ax  = axes[oi]
    sub = meta[meta['Location'] == loc]

    x_pos = 0
    xtick_pos, xtick_lab, xtick_ab = [], [], []
    bar_w = 0.72

    for ab in AB_ORDER:
        for cond in ['HEALTHY', 'DYSBIOTIC']:
            sids = sub[(sub['Antibody'] == ab) & (sub['Condition'] == cond)].index
            sids = [s for s in sids if s in DATA.index]
            if not sids:
                continue

            mean_ra = DATA.loc[sids].mean(axis=0)
            bottom  = 0.0
            for taxon in taxa_list:
                val = mean_ra.get(taxon, 0)
                if val > 0:
                    ax.bar(x_pos, val * 100, width=bar_w, bottom=bottom * 100,
                           color=taxa_pal[taxon], edgecolor='white', linewidth=0.3)
                    bottom += val

            # Store n value in x-axis label
            cond_name = 'Healthy' if cond == 'HEALTHY' else 'Dysbiotic'
            n_sids = len(sids)

            xtick_pos.append(x_pos)
            xtick_lab.append(f"{cond_name}\n(n={n_sids})")
            xtick_ab.append(ab)
            x_pos += 1
        x_pos += 0.4  # gap between fractions

    ax.set_xticks(xtick_pos)
    ax.set_xticklabels(xtick_lab, fontsize=8)

    # Apply colors to xticklabels based on antibody
    for label, ab in zip(ax.get_xticklabels(), xtick_ab):
        label.set_color(AB_C[ab])

    ax.set_ylim(0, 100)
    ax.set_title(loc.title(), fontsize=13, fontweight='bold',
                 color=LOC_C[loc], pad=8)
    ax.set_ylabel('Mean relative abundance (%)', fontsize=11)
    ax.spines[['top', 'right']].set_visible(False)
    ax.yaxis.grid(True, alpha=0.2, linestyle='--')
    ax.set_axisbelow(True)

    # Fraction group brackets
    tick_arr = np.array(xtick_pos)
    ab_lab_done = set()
    for ab in AB_ORDER:
        # Get all indices for this antibody
        ab_idxs = [i for i, a in enumerate(xtick_ab) if a == ab]
        if ab_idxs and ab not in ab_lab_done:
            xm = (tick_arr[ab_idxs[0]] + tick_arr[ab_idxs[-1]]) / 2
            ax.text(xm, -5, AB_LABELS[ab], ha='center', va='top', fontsize=9,
                    fontweight='bold', color=AB_C[ab])
            ab_lab_done.add(ab)

# Legend — top 20 most abundant taxa
top20 = DATA.mean(axis=0).sort_values(ascending=False).head(20).index.tolist()
handles = [mpatches.Patch(color=taxa_pal[t], label=t) for t in top20]
fig.legend(handles=handles, loc='lower center', ncol=10, fontsize=8,
           frameon=True, title='Bacterial taxon (top 20 shown)',
           title_fontsize=9, bbox_to_anchor=(0.5, -0.07))

plt.tight_layout()
fig.savefig(os.path.join(FIG_DIR, 'SupplementaryFigure13.png'),
            dpi=600, bbox_inches='tight', facecolor='white')
plt.close()
print("Saved SupplementaryFigure13.png")
