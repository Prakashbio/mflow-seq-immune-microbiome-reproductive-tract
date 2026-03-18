import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import FuncFormatter
from matplotlib.patches import Ellipse
from scipy.stats import mannwhitneyu
import warnings; warnings.filterwarnings('ignore')

from shared_utils import (load_data, shannon, hellinger, bray_curtis, pcoa,
                           permanova, pct_label,
                           LOC_C, COND_C, LOC_ORDER, FIG_DIR)

otu_rel, meta, tax = load_data()
DATA  = otu_rel.T
taxa  = otu_rel.index.tolist()
PSEUDO = 1e-5

def plot_ellipse(ax, pts, color, n_std=1.96, alpha=0.15):
    if len(pts) < 3: return
    mean = pts.mean(axis=0); cov = np.cov(pts.T)
    if cov.ndim < 2 or np.any(np.isnan(cov)): return
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]; vals = vals[order]; vecs = vecs[:, order]
    angle = np.degrees(np.arctan2(*vecs[:, 0][::-1]))
    w, h = 2 * n_std * np.sqrt(np.abs(vals))
    ax.add_patch(Ellipse(xy=mean, width=w, height=h, angle=angle,
                         facecolor=color, alpha=alpha, edgecolor=color, lw=1.5))

# ── Per-sample diversity ──────────────────────────────────────────────────────
div_rows = [{'SampleID': sid, 'Shannon': shannon(DATA.loc[sid].values)}
            for sid in DATA.index]
div = pd.DataFrame(div_rows).set_index('SampleID').join(meta)

KEY_TAXA = ['Lactobacillus iners', 'Lactobacillus crispatus',
            'Gardnerella vaginalis', 'Sneathia vaginalis',
            'Fannyhessea vaginae', 'Prevotella buccalis']


def significance_label(p):
    label = pct_label(p)
    return 'n.s.' if label == 'ns' else label

# ─────────────────────────────────────────────────────────────────────────────
# MAIN FIGURE
# ─────────────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(22, 20))
fig.patch.set_facecolor('white')


# Uniform vertical spacing across rows.
gs = fig.add_gridspec(3, 3, hspace=0.20, wspace=0.15)
a_axes = []
b_axes = []
c_axes = []
a_cond_handles = None

# ── Panels A: Key taxa abundance Healthy vs Dysbiotic per organ ───────────────
for oi, loc in enumerate(LOC_ORDER):
    ax   = fig.add_subplot(gs[0, oi])
    a_axes.append(ax)
    sub  = meta[meta['Location'] == loc]
    # Use pre-sorted fraction for cleanest comparison
    pre  = sub[sub['Antibody'] == 'PRE']

    x_pos  = 0
    xtick_pos, xtick_lab = [], []
    spp_centers, spp_labels, spp_colors = [], [], []
    panel_max = 0.0
    panel_sig_max = 0.0
    spp_palette = plt.cm.tab10(np.linspace(0, 1, max(1, len(KEY_TAXA))))

    for ti, taxon in enumerate(KEY_TAXA):
        if taxon not in DATA.columns: continue
        pair_positions = {}
        pair_max = 0.0
        for ci, cond in enumerate(['HEALTHY', 'DYSBIOTIC']):
            sids = pre[pre['Condition'] == cond].index
            sids = [s for s in sids if s in DATA.index]
            if not sids: continue
            vals = DATA.loc[sids, taxon].values * 100
            pair_positions[cond] = x_pos
            pair_max = max(pair_max, float(np.nanmax(vals)) if len(vals) else 0.0)
            col  = COND_C[cond]
            bp   = ax.boxplot(vals, positions=[x_pos], widths=0.38,
                              patch_artist=True, showfliers=False,
                              medianprops={'color': 'white', 'lw': 2},
                              whiskerprops={'color': col, 'lw': 1.5},
                              capprops={'color': col, 'lw': 1.5},
                              boxprops={'facecolor': col, 'alpha': 0.50,
                                        'edgecolor': col, 'lw': 1.2})
            np.random.seed(ti * 5 + ci)
            ax.scatter(x_pos + np.random.normal(0, 0.04, len(vals)),
                       vals, c=col, s=35, alpha=0.8, zorder=8,
                       edgecolors='white', lw=0.3)

            xtick_pos.append(x_pos)
            xtick_lab.append('H' if cond == 'HEALTHY' else 'D')
            x_pos += 0.5

        if 'HEALTHY' in pair_positions and 'DYSBIOTIC' in pair_positions:
            vals_h = DATA.loc[pre[pre['Condition']=='HEALTHY'].index.intersection(DATA.index), taxon].values * 100
            vals_d = DATA.loc[pre[pre['Condition']=='DYSBIOTIC'].index.intersection(DATA.index), taxon].values * 100
            if len(vals_h) > 1 and len(vals_d) > 1:
                try:
                    _, p = mannwhitneyu(vals_h, vals_d, alternative='two-sided')
                    y_base = max(pair_max, 0.2)
                    y_bar = y_base * 1.30
                    y_txt = y_base * 1.48
                    hx = pair_positions['HEALTHY']
                    dx = pair_positions['DYSBIOTIC']
                    ax.plot([hx, dx], [y_bar, y_bar], color='#555', lw=1.2, zorder=9)
                    ax.text((hx + dx) / 2, y_txt, significance_label(p),
                            ha='center', va='bottom', fontsize=9, color='#333', zorder=10)
                    panel_sig_max = max(panel_sig_max, y_txt)
                except:
                    pass

        if pair_positions:
            pos_vals = list(pair_positions.values())
            spp_centers.append(float(np.mean(pos_vals)))
            spp_labels.append(' '.join(taxon.split(' ')[1:]))
            spp_colors.append(spp_palette[ti % len(spp_palette)])
            panel_max = max(panel_max, pair_max)

        x_pos += 1

    ax.set_xticks(xtick_pos)
    ax.set_xticklabels(xtick_lab, fontsize=9, fontweight='bold')
    if spp_centers:
        secax = ax.secondary_xaxis('bottom')
        secax.set_xticks(spp_centers)
        secax.set_xticklabels(spp_labels, fontsize=8, fontstyle='italic', fontweight='bold')
        secax.tick_params(axis='x', length=0, pad=18)
        secax.spines['bottom'].set_visible(False)
        for lbl, s_color in zip(secax.get_xticklabels(), spp_colors):
            lbl.set_color(s_color)
    if xtick_pos:
        ax.set_xlim(min(xtick_pos) - 0.6, max(xtick_pos) + 0.6)
    ax.set_ylabel('Relative abundance (%)', fontsize=10)
    ax.set_title(f'A   {loc.title()} — key taxa', fontsize=11,
                 fontweight='bold', color=LOC_C[loc], loc='left', pad=6)
    ax.spines[['top', 'right']].set_visible(False)
    ax.yaxis.grid(True, alpha=0.2, linestyle='--')
    ax.set_axisbelow(True)
    ax.set_yscale('symlog', linthresh=1)
    ax.set_ylim(bottom=0)
    if panel_max > 0:
        y_top = max(panel_max * 1.7, panel_sig_max * 1.15, ax.get_ylim()[1])
        ax.set_ylim(top=y_top)
    y_max = ax.get_ylim()[1]
    yticks = [0, 1]
    p = 2
    while 10 ** p <= y_max * 1.02:
        yticks.append(10 ** p)
        p += 1
    ax.set_yticks(yticks)
    ax.yaxis.set_major_formatter(FuncFormatter(
        lambda y, _: '0' if np.isclose(y, 0) else (f"{int(y)}" if y >= 1 else '')
    ))

if a_axes:
    a_cond_handles = [
        mpatches.Patch(facecolor=COND_C['HEALTHY'], edgecolor=COND_C['HEALTHY'], alpha=0.55, label='Healthy'),
        mpatches.Patch(facecolor=COND_C['DYSBIOTIC'], edgecolor=COND_C['DYSBIOTIC'], alpha=0.55, label='Dysbiotic')
    ]

# ── Panels B: Shannon diversity Healthy vs Dysbiotic per organ ───────────────
for oi, loc in enumerate(LOC_ORDER):
    ax  = fig.add_subplot(gs[1, oi])
    b_axes.append(ax)
    sub = div[(div['Location'] == loc) & (div['Antibody'] == 'PRE')]
    counts = {}

    for ci, cond in enumerate(['HEALTHY', 'DYSBIOTIC']):
        vals = sub[sub['Condition'] == cond]['Shannon'].dropna()
        counts[cond] = len(vals)
        if len(vals) == 0: continue
        col  = COND_C[cond]
        ax.boxplot(vals, positions=[ci], widths=0.38, patch_artist=True,
                   showfliers=False,
                   medianprops={'color': 'white', 'lw': 2.5},
                   whiskerprops={'color': col, 'lw': 1.5},
                   capprops={'color': col, 'lw': 1.5},
                   boxprops={'facecolor': col, 'alpha': 0.55,
                              'edgecolor': col, 'lw': 1.2})
        np.random.seed(ci * 100)
        ax.scatter(ci + np.random.normal(0, 0.06, len(vals)), vals,
                   c=col, s=50, alpha=0.85, zorder=8,
                   edgecolors='white', lw=0.4)

    # MWU
    h = sub[sub['Condition']=='HEALTHY']['Shannon'].dropna()
    d = sub[sub['Condition']=='DYSBIOTIC']['Shannon'].dropna()
    if len(h) > 1 and len(d) > 1:
        try:
            _, p = mannwhitneyu(h, d, alternative='two-sided')
            ymax = max(h.max(), d.max())
            spread = max(float(sub['Shannon'].max() - sub['Shannon'].min()), 0.25)
            y_bar = ymax + 0.14 * spread
            y_txt = y_bar + 0.08 * spread
            ax.annotate('', xy=(1, y_bar), xytext=(0, y_bar),
                        arrowprops=dict(arrowstyle='-', color='#555', lw=1.5))
            ax.text(0.5, y_txt, significance_label(p), ha='center', fontsize=11)
            ax.set_ylim(top=max(ax.get_ylim()[1], y_txt + 0.15 * spread))
        except: pass

    ax.set_xticks([0, 1])
    ax.set_xticklabels([
        f"Healthy\n(n={counts.get('HEALTHY', 0)})",
        f"Dysbiotic\n(n={counts.get('DYSBIOTIC', 0)})"
    ], fontsize=11, fontweight='bold')
    ax.get_xticklabels()[0].set_color(COND_C['HEALTHY'])
    ax.get_xticklabels()[1].set_color(COND_C['DYSBIOTIC'])
    ax.set_ylabel("Shannon H'", fontsize=11)
    ax.set_title(f'B   {loc.title()} — Shannon diversity', fontsize=11,
                 fontweight='bold', color=LOC_C[loc], loc='left', pad=6)
    ax.spines[['top', 'right']].set_visible(False)
    ax.yaxis.grid(True, alpha=0.22, linestyle='--')
    ax.set_axisbelow(True)
    ax.set_ylim(bottom=0)

# ── Panels C: PCoA per organ coloured by condition ───────────────────────────
for oi, loc in enumerate(LOC_ORDER):
    ax   = fig.add_subplot(gs[2, oi])
    c_axes.append(ax)
    sids = meta[(meta['Location'] == loc) & (meta['Antibody'] == 'PRE')].index
    sids = [s for s in sids if s in DATA.index]
    if len(sids) < 4:
        ax.text(0.5, 0.5, 'Insufficient data', transform=ax.transAxes, ha='center')
        continue

    Xi   = otu_rel[sids].T.values
    Xh   = hellinger(Xi)
    Di   = bray_curtis(Xh)
    ci_c, vi = pcoa(Di, n_axes=2)
    conds = meta.loc[sids, 'Condition'].values
    legend_handles = []

    for cond in ['HEALTHY', 'DYSBIOTIC']:
        idx = np.where(conds == cond)[0]
        if len(idx) == 0:
            continue
        pts = ci_c[idx, :2]
        col = COND_C[cond]
        scatter = ax.scatter(pts[:, 0], pts[:, 1], c=col, s=65, alpha=0.85,
                             edgecolors='white', lw=0.5,
                             label=f"{cond.title()} (n={len(idx)})", zorder=5)
        legend_handles.append(scatter)
        plot_ellipse(ax, pts, col)

    pr = permanova(Di, conds, n_perm=499)
    ax.text(0.97, 0.03,
            f"R²={pr['R2']:.3f}, p={pr['p']:.3f} {pct_label(pr['p'])}",
            transform=ax.transAxes, ha='right', va='bottom', fontsize=9,
            style='italic', color='#333',
            bbox=dict(boxstyle='round,pad=0.3', fc='#f8f8f8', ec='#ccc', alpha=0.85))
    ax.set_title(f'C   {loc.title()} — PCoA (Pre-sorted)', fontsize=11,
                 fontweight='bold', color=LOC_C[loc], loc='left', pad=6)
    ax.set_xlabel(f'PCoA1 ({vi[0]*100:.1f}%)', fontsize=9)
    ax.set_ylabel(f'PCoA2 ({vi[1]*100:.1f}%)', fontsize=9)
    ax.spines[['top', 'right']].set_visible(False)
    ax.axhline(0, color='#ddd', lw=0.8, linestyle='--')
    ax.axvline(0, color='#ddd', lw=0.8, linestyle='--')
    if legend_handles:
        ax.legend(handles=legend_handles, loc='upper center', bbox_to_anchor=(0.5, -0.07),
                  ncol=max(1, len(legend_handles)), fontsize=9, frameon=False,
                  columnspacing=1.6, handletextpad=0.6)

plt.tight_layout(rect=[0.02, 0.07, 0.98, 0.98])

# Dedicated legend axis under middle Panel A for robust placement.
if a_cond_handles and a_axes and b_axes:
    a_mid = a_axes[1] if len(a_axes) > 1 else a_axes[0]
    b_mid = b_axes[1] if len(b_axes) > 1 else b_axes[0]
    a_pos = a_mid.get_position()
    b_pos = b_mid.get_position()

    leg_h = 0.03
    leg_y = max(b_pos.y1 + 0.002, (a_pos.y0 + b_pos.y1) / 2 - (leg_h / 2))
    legend_ax = fig.add_axes([a_pos.x0, leg_y, a_pos.width, leg_h])
    legend_ax.axis('off')
    legend_ax.legend(handles=a_cond_handles, loc='center', ncol=2, frameon=False,
                     fontsize=10, columnspacing=1.8, handlelength=1.4)

fig.savefig(os.path.join(FIG_DIR, 'Figure6.png'),
            dpi=600, bbox_inches='tight', facecolor='white')
plt.close()
print("Saved Figure6.png")
