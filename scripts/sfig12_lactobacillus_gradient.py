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

from shared_utils import (load_data, pct_label,
                           LOC_C, AB_C, COND_C, LOC_ORDER, AB_ORDER, AB_LABELS, FIG_DIR)

otu_rel, meta, tax = load_data()
DATA   = otu_rel.T
taxa   = otu_rel.index.tolist()
PSEUDO = 1e-5

lacto_taxa   = [t for t in taxa if 'Lactobacillus' in t]
dysbio_taxa  = ['Gardnerella vaginalis', 'Gardnerella leopoldii',
                'Sneathia vaginalis', 'Fannyhessea vaginae',
                'Prevotella buccalis', 'Prevotella oris',
                'Parvimonas micra', 'Metamycoplasma hominis']
dysbio_taxa  = [t for t in dysbio_taxa if t in taxa]

LACTO_COLORS = {
    'Lactobacillus crispatus':  '#1565C0',
    'Lactobacillus iners':      '#2196F3',
    'Lactobacillus gasseri':    '#64B5F6',
    'Lactobacillus jensenii':   '#90CAF9',
    'Lactobacillus delbrueckii':'#BBDEFB',
}

LACTO_HATCHES = {
    'Lactobacillus crispatus': '',
    'Lactobacillus iners': '///',
    'Lactobacillus gasseri': '\\\\\\',
    'Lactobacillus jensenii': 'xxx',
    'Lactobacillus delbrueckii': '...',
}

def group_mean(loc, ab, condition=None):
    q = meta[(meta['Location']==loc) & (meta['Antibody']==ab)]
    if condition: q = q[q['Condition']==condition]
    sids = [s for s in q.index if s in DATA.index]
    if not sids: return pd.Series(0.0, index=taxa)
    return DATA.loc[sids].mean(axis=0)

# ─────────────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(22, 20))
fig.patch.set_facecolor('white')

gs = fig.add_gridspec(3, 3, hspace=0.25, wspace=0.15)

# ── Panel A: Lactobacillus species stacked bars ───────────────────────────────
ax_a = fig.add_subplot(gs[0, :])

x_pos = 0
xtick_pos, xtick_lab = [], []
bar_w = 0.72

for loc in LOC_ORDER:
    for ab in AB_ORDER:
        mu = group_mean(loc, ab)
        bottom = 0.0
        for lt in lacto_taxa:
            val = mu.get(lt, 0) * 100
            if val > 0:
                ax_a.bar(x_pos, val, width=bar_w, bottom=bottom,
                         color=LACTO_COLORS.get(lt, '#888888'),
                         edgecolor='white', linewidth=0.4)
                if val > 5:
                    ax_a.text(x_pos, bottom + val/2, lt.split(' ')[1][:6],
                              ha='center', va='center', fontsize=6.5,
                              color='white', fontweight='bold', style='italic')
                bottom += val
        # Non-lacto total as light grey
        non_lacto = 100 - bottom
        if non_lacto > 0:
            ax_a.bar(x_pos, non_lacto, width=bar_w, bottom=bottom,
                     color='#DDDDDD', edgecolor='white', linewidth=0.4)
        n_sids = len([s for s in meta[(meta['Location']==loc)&(meta['Antibody']==ab)].index
                      if s in DATA.index])
        xtick_pos.append(x_pos)
        xtick_lab.append(f'{AB_LABELS[ab]}\n(n={n_sids})')
        x_pos += 1
    x_pos += 0.5

ax_a.set_xticks(xtick_pos)
ax_a.set_xticklabels(xtick_lab, fontsize=9, fontweight='bold')
ax_a.set_ylim(0, 115)
ax_a.set_ylabel('Mean relative abundance (%)', fontsize=11)
ax_a.set_title('A   Lactobacillus species composition by anatomical site and Ig fraction (grey = non-Lactobacillus taxa)',
               fontsize=11, fontweight='bold', loc='left', pad=6)
ax_a.spines[['top', 'right']].set_visible(False)
ax_a.yaxis.grid(True, alpha=0.22, linestyle='--')
ax_a.set_axisbelow(True)

# Add organ labels
for oi, loc in enumerate(LOC_ORDER):
    xm = oi * (len(AB_ORDER) + 0.5) + (len(AB_ORDER) - 1) / 2
    ax_a.text(xm, 105, loc.title(), ha='center', va='bottom',
              fontsize=12, fontweight='bold', color=LOC_C[loc])

# Legend
handles_a = [mpatches.Patch(color=LACTO_COLORS.get(lt, '#888888'), label=lt)
             for lt in lacto_taxa]
handles_a += [mpatches.Patch(color='#DDDDDD', label='Non-Lactobacillus')]
ax_a.legend(handles=handles_a, fontsize=8.5, frameon=True, ncol=6,
            loc='lower center', bbox_to_anchor=(0.5, -0.16))

# ── Panels B: Log2FC vs pre-sorted per Lactobacillus species ─────────────────
for oi, loc in enumerate(LOC_ORDER):
    ax_b = fig.add_subplot(gs[1, oi])
    mu_pre = group_mean(loc, 'PRE')

    for ai, ab in enumerate(['A', 'M', 'G']):
        mu_ab  = group_mean(loc, ab)
        lfc    = {lt: np.log2((mu_ab.get(lt, 0) + PSEUDO) /
                               (mu_pre.get(lt, 0) + PSEUDO))
                  for lt in lacto_taxa}
        xbase  = ai * (len(lacto_taxa) + 1.5)
        ab_col = {'A': '#FFB366', 'M': '#B39DDB', 'G': '#81C784'}[ab]  # light colors
        for ti, lt in enumerate(lacto_taxa):
            val = lfc[lt]
            hatch = LACTO_HATCHES.get(lt, '')
            ax_b.bar(xbase + ti, val, width=0.75, color=ab_col, alpha=0.8,
                     edgecolor='white', linewidth=0.4, hatch=hatch)

    ax_b.axhline(0, color='#333', lw=1.5, linestyle='--', alpha=0.7)
    n_ab = 3
    n_lt = len(lacto_taxa)
    tick_pos_b = [(i * (n_lt + 1.5) + (n_lt - 1) / 2) for i in range(n_ab)]
    ax_b.set_xticks(tick_pos_b)
    ax_b.set_xticklabels(['IgA', 'IgM', 'IgG'], fontsize=10, fontweight='bold')
    ax_b.get_xticklabels()[0].set_color(AB_C['A'])
    ax_b.get_xticklabels()[1].set_color(AB_C['M'])
    ax_b.get_xticklabels()[2].set_color(AB_C['G'])
    ax_b.set_ylabel('log₂(fraction / pre-sorted)', fontsize=10)
    ax_b.set_title(f'B   {loc.title()} — Lactobacillus enrichment vs pre-sorted',
                   fontsize=11, fontweight='bold', color=LOC_C[loc], loc='left', pad=6)
    ax_b.spines[['top', 'right']].set_visible(False)
    ax_b.yaxis.grid(True, alpha=0.22, linestyle='--')
    ax_b.set_axisbelow(True)
    # Legend: antibody colors + species patterns (horizontal center below graph) - only on first panel
    if oi == 1:
        handles_b_ab = [mpatches.Patch(color='#FFB366', alpha=0.8, label='IgA'),
                        mpatches.Patch(color='#B39DDB', alpha=0.8, label='IgM'),
                        mpatches.Patch(color='#81C784', alpha=0.8, label='IgG')]
        handles_b_sp = [mpatches.Patch(facecolor='white', edgecolor='black', hatch=LACTO_HATCHES.get(lt, ''),
                                       label=lt, linewidth=0.5)
                       for lt in lacto_taxa]
        all_handles_b = handles_b_ab + handles_b_sp
        ax_b.legend(handles=all_handles_b, fontsize=8, frameon=True, ncol=8,
                    loc='lower center', bbox_to_anchor=(0.5, -0.15))

# ── Panel C: Dysbiotic taxa coating by condition ──────────────────────────────
for oi, loc in enumerate(LOC_ORDER):
    ax_c = fig.add_subplot(gs[2, oi])

    x_p = 0
    ytick_p, ytick_l = [], []
    for ti, taxon in enumerate(dysbio_taxa):
        if taxon not in taxa: continue
        for ci, cond in enumerate(['HEALTHY', 'DYSBIOTIC']):
            a_sids = [s for s in meta[(meta['Location']==loc)&(meta['Antibody']=='A')&
                                       (meta['Condition']==cond)].index if s in DATA.index]
            m_sids = [s for s in meta[(meta['Location']==loc)&(meta['Antibody']=='M')&
                                       (meta['Condition']==cond)].index if s in DATA.index]
            if not a_sids or not m_sids:
                x_p += 1; continue
            mu_a = DATA.loc[a_sids, taxon].mean() * 100
            mu_m = DATA.loc[m_sids, taxon].mean() * 100
            col  = COND_C[cond]
            # bar pair: IgA and IgM side by side
            ax_c.barh(x_p,        mu_a, height=0.35, color='#E66101', alpha=0.7 if cond=='HEALTHY' else 0.95, edgecolor='white')
            ax_c.barh(x_p - 0.35, mu_m, height=0.35, color='#5E3C99', alpha=0.7 if cond=='HEALTHY' else 0.95, edgecolor='white', hatch='///' if cond=='DYSBIOTIC' else '')
            ytick_p.append(x_p - 0.175)
            ytick_l.append(f"{taxon.split(' ')[1][:6]} ({'H' if cond=='HEALTHY' else 'D'})")
            x_p += 1
        x_p += 0.3

    ax_c.set_yticks(ytick_p)
    ax_c.set_yticklabels(ytick_l, fontsize=7.5, style='italic')
    ax_c.invert_yaxis()
    ax_c.set_xlabel('Mean rel. abundance (%)', fontsize=10)
    ax_c.set_title(f'C   {loc.title()} — dysbiotic taxa IgA/IgM by condition',
                   fontsize=11, fontweight='bold', color=LOC_C[loc], loc='left', pad=6)
    ax_c.spines[['top', 'right']].set_visible(False)
    ax_c.xaxis.grid(True, alpha=0.22, linestyle='--')
    ax_c.set_axisbelow(True)
    if oi == 1:
        handles_c = [mpatches.Patch(color='#E66101', alpha=0.8, label='IgA fraction'),
                     mpatches.Patch(color='#5E3C99', alpha=0.8, label='IgM fraction'),
                     mpatches.Patch(color='#888', alpha=0.5, label='Healthy'),
                     mpatches.Patch(color='#888', alpha=1.0, hatch='///', label='Dysbiotic')]
        ax_c.legend(handles=handles_c, fontsize=8, frameon=True, ncol=4,
                    loc='lower center', bbox_to_anchor=(0.5, -0.18))

plt.tight_layout()
fig.savefig(os.path.join(FIG_DIR, 'SupplementaryFigure12.png'),
            dpi=600, bbox_inches='tight', facecolor='white')
plt.close()
print("Saved SupplementaryFigure12.png")
