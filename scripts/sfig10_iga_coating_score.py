import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import mannwhitneyu
import warnings; warnings.filterwarnings('ignore')

from shared_utils import load_data, pct_label, LOC_C, COND_C, LOC_ORDER

NEW_FIG_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'new_figures')
os.makedirs(NEW_FIG_DIR, exist_ok=True)

otu_rel, meta, tax = load_data()
DATA   = otu_rel.T
taxa   = otu_rel.index.tolist()
PSEUDO = 1e-6

# Custom diverging colourmap: purple (IgM) → white → orange (IgA)
iga_cmap = LinearSegmentedColormap.from_list(
    'iga_coating', ['#5E3C99', '#F7F7F7', '#E66101'], N=256)

# ── Compute per-taxon IgA coating score ──────────────────────────────────────
def iga_coating_score(loc, ab_a='A', ab_m='M'):
    """IgA / (IgA + IgM) mean abundance per taxon for a given location."""
    a_sids = [s for s in meta[(meta['Location']==loc)&(meta['Antibody']==ab_a)].index
              if s in DATA.index]
    m_sids = [s for s in meta[(meta['Location']==loc)&(meta['Antibody']==ab_m)].index
              if s in DATA.index]
    mu_a = DATA.loc[a_sids].mean(axis=0) if a_sids else pd.Series(0, index=taxa)
    mu_m = DATA.loc[m_sids].mean(axis=0) if m_sids else pd.Series(0, index=taxa)
    score = pd.Series(np.nan, index=taxa)
    for t in taxa:
        tot = mu_a[t] + mu_m[t]
        if tot > PSEUDO:
            score[t] = mu_a[t] / tot
    return score, mu_a, mu_m

# ─────────────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(22, 18))
fig.patch.set_facecolor('white')

gs_a  = fig.add_gridspec(1, 1, left=0.15, right=0.87, top=0.95, bottom=0.42)
gs_bc = fig.add_gridspec(1, 4, left=0.06, right=0.96, top=0.38, bottom=0.06,
                         width_ratios=[1.2, 1, 1, 1], wspace=0.40)

# ── Panel A: dot-plot heatmap of IgA coating score ───────────────────────────
ax_a = fig.add_subplot(gs_a[0, 0])

score_mat = np.full((len(taxa), len(LOC_ORDER)), np.nan)
size_mat  = np.zeros((len(taxa), len(LOC_ORDER)))
for oi, loc in enumerate(LOC_ORDER):
    sc, mu_a, mu_m = iga_coating_score(loc)
    for ti, t in enumerate(taxa):
        score_mat[ti, oi] = sc[t]
        size_mat[ti, oi]  = (mu_a[t] + mu_m[t]) * 100  # mean %

# Sort taxa by mean IgA coating score descending
mean_score = np.nanmean(score_mat, axis=1)
sort_idx   = np.argsort(mean_score)[::-1]
score_mat  = score_mat[sort_idx]
size_mat   = size_mat[sort_idx]
taxa_sorted = [taxa[i] for i in sort_idx]

x_positions = np.arange(len(LOC_ORDER)) * 3
for ti, taxon in enumerate(taxa_sorted):
    for oi, loc in enumerate(LOC_ORDER):
        s  = score_mat[ti, oi]
        sz = size_mat[ti, oi]
        if np.isnan(s): continue
        color = iga_cmap(s)
        circle_size = np.clip(sz * 60 + 20, 20, 500)
        ax_a.scatter(x_positions[oi], ti, s=circle_size, c=[color],
                     edgecolors='#444', linewidths=0.5, zorder=5)
        if circle_size > 80:
            ax_a.text(x_positions[oi], ti, f'{s:.2f}', ha='center', va='center',
                      fontsize=6, color='white' if (s < 0.25 or s > 0.75) else '#333',
                      fontweight='bold')

ax_a.set_xticks(x_positions)
ax_a.set_xticklabels([l.title() for l in LOC_ORDER], fontsize=12, fontweight='bold')
for tick, loc in zip(ax_a.get_xticklabels(), LOC_ORDER):
    tick.set_color(LOC_C[loc])
ax_a.set_yticks(range(len(taxa_sorted)))
ax_a.set_yticklabels(taxa_sorted, fontsize=8.5, style='italic')
ax_a.set_xlim(-2, x_positions[-1] + 2)
ax_a.set_ylim(-1, len(taxa_sorted))
ax_a.set_title('A   IgA coating score per taxon per anatomical site  '
               '(circle size = mean relative abundance; colour = score)',
               fontsize=11, fontweight='bold', loc='left', pad=6)
ax_a.spines[['top', 'right', 'bottom', 'left']].set_visible(False)
ax_a.xaxis.grid(True, alpha=0.15, linestyle='--')
ax_a.set_axisbelow(True)

# Colourbar
sm = plt.cm.ScalarMappable(cmap=iga_cmap, norm=plt.Normalize(0, 1))
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax_a, orientation='vertical', fraction=0.015, pad=0.01)
cbar.set_label('IgA coating score [IgA/(IgA+IgM)]', fontsize=9)
cbar.set_ticks([0, 0.33, 0.5, 0.67, 1])
cbar.set_ticklabels(['IgM\ndominant', '0.33', '0.50\n(equal)', '0.67', 'IgA\ndominant'])

# ── Panel B: Lactobacillus vs other genera ────────────────────────────────────
ax_b = fig.add_subplot(gs_bc[0, 0])

lacto_taxa = [t for t in taxa if 'Lactobacillus' in t or 'Lacticasibacillus' in t]
other_taxa = [t for t in taxa if t not in lacto_taxa]

for oi, loc in enumerate(LOC_ORDER):
    sc, _, _ = iga_coating_score(loc)
    lacto_scores = sc[lacto_taxa].dropna().values
    other_scores = sc[other_taxa].dropna().values

    xbase = oi * 2.5
    for xi, (vals, label, col) in enumerate([
        (lacto_scores, 'Lactobacillus', '#4DAC26'),
        (other_scores, 'Other genera',  '#D7191C')
    ]):
        if len(vals) == 0: continue
        bp = ax_b.boxplot(vals, positions=[xbase + xi * 1.0], widths=0.75,
                          patch_artist=True, showfliers=True,
                          medianprops={'color': 'white', 'lw': 2.5},
                          whiskerprops={'color': col, 'lw': 1.5},
                          capprops={'color': col, 'lw': 1.5},
                          flierprops={'marker': 'o', 'ms': 4, 'alpha': 0.5,
                                      'markerfacecolor': col},
                          boxprops={'facecolor': col, 'alpha': 0.50,
                                     'edgecolor': col, 'lw': 1.2})
        np.random.seed(oi * 10 + xi)
        ax_b.scatter(xbase + xi * 1.0 + np.random.normal(0, 0.07, len(vals)),
                     vals, c=col, s=30, alpha=0.8, zorder=8,
                     edgecolors='white', lw=0.3)

    if len(lacto_scores) > 1 and len(other_scores) > 1:
        try:
            _, p = mannwhitneyu(lacto_scores, other_scores, alternative='two-sided')
            ymax = max(lacto_scores.max() if len(lacto_scores) else 0,
                       other_scores.max() if len(other_scores) else 0)
            ax_b.text(xbase + 0.5, min(ymax + 0.06, 1.05),
                      pct_label(p), ha='center', fontsize=11, fontweight='bold')
        except: pass

ax_b.axhline(0.5, color='#888', lw=1.2, linestyle=':', alpha=0.7, label='Equal coating (0.5)')
ax_b.set_xticks([0.5, 3.0, 5.5])
ax_b.set_xticklabels([l.title() for l in LOC_ORDER], fontsize=10, fontweight='bold')
for tick, loc in zip(ax_b.get_xticklabels(), LOC_ORDER):
    tick.set_color(LOC_C[loc])
ax_b.set_ylabel('IgA coating score [IgA/(IgA+IgM)]', fontsize=10)
ax_b.set_title('B   IgA coating: Lactobacillus vs other genera',
               fontsize=11, fontweight='bold', loc='left', pad=6)
ax_b.spines[['top', 'right']].set_visible(False)
ax_b.yaxis.grid(True, alpha=0.22, linestyle='--')
ax_b.set_axisbelow(True)
ax_b.set_ylim(-0.05, 1.1)
handles_b = [mpatches.Patch(color='#4DAC26', alpha=0.6, label='Lactobacillus spp.'),
             mpatches.Patch(color='#D7191C', alpha=0.6, label='Other genera')]
ax_b.legend(handles=handles_b, fontsize=9, frameon=True,
            loc='lower center', bbox_to_anchor=(0.5, -0.1), ncols=2)

# ── Panel C: IgA coating score by condition per organ ────────────────────────
for oi, loc in enumerate(LOC_ORDER):
    ax_c = fig.add_subplot(gs_bc[0, oi + 1])
    # per-sample IgA coating score: sum(IgA) / (sum(IgA) + sum(IgM)) per sample
    patient_scores = {'HEALTHY': [], 'DYSBIOTIC': []}
    meta_loc = meta[meta['Location'] == loc]
    for sid_a in meta_loc[meta_loc['Antibody'] == 'A'].index:
        if sid_a not in DATA.index: continue
        pat  = sid_a.split('_')[0]
        cond = meta.loc[sid_a, 'Condition']
        # Find matching IgM sample for same patient
        m_match = meta_loc[(meta_loc['Antibody'] == 'M') &
                            (meta_loc.index.str.startswith(pat + '_'))].index
        if not len(m_match): continue
        sid_m = m_match[0]
        if sid_m not in DATA.index: continue
        a_tot = DATA.loc[sid_a].sum()
        m_tot = DATA.loc[sid_m].sum()
        if a_tot + m_tot > 0:
            patient_scores[cond].append(a_tot / (a_tot + m_tot))

    for ci, cond in enumerate(['HEALTHY', 'DYSBIOTIC']):
        vals = patient_scores[cond]
        if not vals: continue
        col  = COND_C[cond]
        ax_c.boxplot(vals, positions=[ci], widths=0.55, patch_artist=True,
                     showfliers=True,
                     medianprops={'color': 'white', 'lw': 2.5},
                     whiskerprops={'color': col, 'lw': 1.5},
                     capprops={'color': col, 'lw': 1.5},
                     flierprops={'marker': 'o', 'ms': 4, 'alpha': 0.6,
                                  'markerfacecolor': col},
                     boxprops={'facecolor': col, 'alpha': 0.55,
                                'edgecolor': col, 'lw': 1.2})
        np.random.seed(ci + oi * 3)
        ax_c.scatter(ci + np.random.normal(0, 0.06, len(vals)), vals,
                     c=col, s=50, alpha=0.85, zorder=8,
                     edgecolors='white', lw=0.4)

    h_v = patient_scores['HEALTHY']
    d_v = patient_scores['DYSBIOTIC']
    if len(h_v) > 1 and len(d_v) > 1:
        try:
            _, p = mannwhitneyu(h_v, d_v, alternative='two-sided')
            ymax = max(max(h_v), max(d_v))
            ax_c.annotate('', xy=(1, ymax + 0.05), xytext=(0, ymax + 0.05),
                          arrowprops=dict(arrowstyle='-', color='#555', lw=1.5))
            ax_c.text(0.5, ymax + 0.07, pct_label(p), ha='center', fontsize=12)
        except: pass

    ax_c.axhline(0.5, color='#888', lw=1.2, linestyle=':', alpha=0.7)
    ax_c.set_xticks([0, 1])
    ax_c.set_xticklabels([
        f'Healthy\n(n={len(patient_scores["HEALTHY"])})',
        f'Dysbiotic\n(n={len(patient_scores["DYSBIOTIC"])})'],
        fontsize=10, fontweight='bold')
    ax_c.get_xticklabels()[0].set_color(COND_C['HEALTHY'])
    ax_c.get_xticklabels()[1].set_color(COND_C['DYSBIOTIC'])
    ax_c.set_ylim(-0.05, 1.15)
    ax_c.set_ylabel('Per-sample IgA coating score', fontsize=9)
    ax_c.set_title(f'C   {loc.title()}', fontsize=11, fontweight='bold',
                   color=LOC_C[loc], pad=12)
    ax_c.spines[['top', 'right']].set_visible(False)
    ax_c.yaxis.grid(True, alpha=0.20, linestyle='--')
    ax_c.set_axisbelow(True)

fig.savefig(os.path.join(NEW_FIG_DIR, 'SupplementaryFigure10.png'),
            dpi=600, bbox_inches='tight', facecolor='white')
plt.close()
print("Saved SupplementaryFigure10.png")
