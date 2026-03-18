import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import warnings; warnings.filterwarnings('ignore')

from shared_utils import load_data, LOC_C, LOC_ORDER, AB_ORDER, AB_LABELS, FIG_DIR

# ── Load ──────────────────────────────────────────────────────────────────────
otu_rel, meta, tax = load_data()
DATA = otu_rel.T  # (samples × taxa)

# ─────────────────────────────────────────────────────────────────────────────
# SUPPLEMENTARY: Phylum-level composition
# ─────────────────────────────────────────────────────────────────────────────
phylum_pal = {
    'Bacillota':       '#4393C3',
    'Actinomycetota':  '#D6604D',
    'Pseudomonadota':  '#74C476',
    'Bacteroidota':    '#FD8D3C',
    'Fusobacteriota':  '#9E9AC8',
    'Mycoplasmatota':  '#FDAE6B',
    'Other':           '#BBBBBB',
}

# Map taxa to phylum
taxon_phylum = tax['Phylum'].to_dict()

fig2, axes2 = plt.subplots(1, 3, figsize=(20, 7))
fig2.patch.set_facecolor('white')

for oi, loc in enumerate(LOC_ORDER):
    ax  = axes2[oi]
    sub = meta[meta['Location'] == loc]
    x_pos = 0
    xtick_pos, xtick_lab = [], []

    for ab in AB_ORDER:
        for cond in ['HEALTHY', 'DYSBIOTIC']:
            sids = sub[(sub['Antibody'] == ab) & (sub['Condition'] == cond)].index
            sids = [s for s in sids if s in DATA.index]
            if not sids:
                continue

            mean_ra = DATA.loc[sids].mean(axis=0)
            phylum_sums = {}
            for taxon, val in mean_ra.items():
                if val > 0:
                    ph = taxon_phylum.get(taxon, 'Other')
                    phylum_sums[ph] = phylum_sums.get(ph, 0) + val

            bottom = 0.0
            for ph, col in phylum_pal.items():
                val = phylum_sums.get(ph, 0)
                if val > 0:
                    ax.bar(x_pos, val * 100, width=0.72, bottom=bottom * 100,
                           color=col, edgecolor='white', linewidth=0.5, label=ph if oi == 0 else '_')
                    if val * 100 > 8:
                        ax.text(x_pos, bottom * 100 + val * 50, ph[:8],
                                ha='center', va='center', fontsize=6.5,
                                color='white', fontweight='bold')
                    bottom += val

            cond_sym = 'H' if cond == 'HEALTHY' else 'D'
            xtick_pos.append(x_pos)
            xtick_lab.append(f"{AB_LABELS[ab]}\n({cond_sym})")
            x_pos += 1
        x_pos += 0.4

    ax.set_xticks(xtick_pos)
    ax.set_xticklabels(xtick_lab, fontsize=8)
    ax.set_ylim(0, 100)
    ax.set_title(loc.title(), fontsize=13, fontweight='bold',
                 color=LOC_C[loc], pad=8)
    ax.set_ylabel('Mean relative abundance (%)', fontsize=11)
    ax.spines[['top', 'right']].set_visible(False)
    ax.yaxis.grid(True, alpha=0.2, linestyle='--')
    ax.set_axisbelow(True)

handles2 = [mpatches.Patch(color=c, label=p) for p, c in phylum_pal.items()]
fig2.legend(handles=handles2, loc='lower center', ncol=10, fontsize=10,
            frameon=True, title='Phylum', title_fontsize=10,
            bbox_to_anchor=(0.5, -0.07))

plt.tight_layout()
fig2.savefig(os.path.join(FIG_DIR, 'SupplementaryFigure1.png'),
             dpi=600, bbox_inches='tight', facecolor='white')
plt.close()
print("Saved SupplementaryFigure1.png")
