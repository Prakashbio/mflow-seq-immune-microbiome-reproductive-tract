import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
from sklearn.decomposition import PCA
import warnings; warnings.filterwarnings('ignore')

from shared_utils import (load_data, clr, bray_curtis, permanova, pct_label,
                           hellinger,
                           LOC_C, AB_C, LOC_ORDER, AB_LABELS, FIG_DIR)

otu_rel, meta, tax = load_data()
DATA  = otu_rel.T
taxa  = otu_rel.index.tolist()

# ─────────────────────────────────────────────────────────────────────────────
# MAIN FIGURE
# ─────────────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(24, 18))
fig.patch.set_facecolor('white')

gs = GridSpec(2, 3, figure=fig, hspace=0.15, wspace=0.40)

# ── Panel A: CLR-PCA biplot ───────────────────────────────────────────────────
ax_pca = fig.add_subplot(gs[0, 0:2])

X_clr = clr(DATA.values, pseudocount=0.5)
pca   = PCA(n_components=4)
scores = pca.fit_transform(X_clr)
loads  = pca.components_   # shape (n_PC, n_taxa)
var_e  = pca.explained_variance_ratio_

# Sample scatter coloured by Location
samples = DATA.index.tolist()
for loc in LOC_ORDER:
    sids = meta[meta['Location'] == loc].index
    idx  = [i for i, s in enumerate(samples) if s in sids]
    pts  = scores[idx, :2]
    ax_pca.scatter(pts[:, 0], pts[:, 1], c=LOC_C[loc], s=55, alpha=0.75,
                   edgecolors='white', lw=0.4, label=loc.title(), zorder=5)

# Loading arrows (top 10 taxa by loading magnitude)
load_mag = np.sqrt(loads[0] ** 2 + loads[1] ** 2)
top_idx  = load_mag.argsort()[::-1][:10]
scale    = np.abs(scores[:, :2]).max() * 0.85
loading_points = []

for i in top_idx:
    lx = loads[0, i] * scale
    ly = loads[1, i] * scale
    ax_pca.annotate('', xy=(lx, ly), xytext=(0, 0),
                    arrowprops=dict(arrowstyle='->', color='#333', lw=1.4,
                                    mutation_scale=12))
    ha = 'left' if lx >= 0 else 'right'
    tx = lx * 1.08
    ty = ly * 1.08
    loading_points.extend([(lx, ly), (tx, ty)])
    ax_pca.text(tx, ty, taxa[i].split(' ')[1],
                fontsize=8, ha=ha, va='center', style='italic', color='#B22222')

if loading_points:
    pca_points = np.vstack([scores[:, :2], np.array(loading_points)])
else:
    pca_points = scores[:, :2]

x_span = max(pca_points[:, 0].max() - pca_points[:, 0].min(), 1e-6)
y_span = max(pca_points[:, 1].max() - pca_points[:, 1].min(), 1e-6)
ax_pca.set_xlim(pca_points[:, 0].min() - 0.08 * x_span,
                pca_points[:, 0].max() + 0.08 * x_span)
ax_pca.set_ylim(pca_points[:, 1].min() - 0.08 * y_span,
                pca_points[:, 1].max() + 0.08 * y_span)

ax_pca.axhline(0, color='#ccc', lw=0.8, linestyle='--')
ax_pca.axvline(0, color='#ccc', lw=0.8, linestyle='--')
ax_pca.set_xlabel(f'PC1 ({var_e[0]*100:.1f}%)', fontsize=11)
ax_pca.set_ylabel(f'PC2 ({var_e[1]*100:.1f}%)', fontsize=11)
ax_pca.set_title('A   CLR-PCA biplot — samples (circles) and species loadings (arrows)',
                 fontsize=11, fontweight='bold', loc='left', pad=6)
ax_pca.spines[['top', 'right']].set_visible(False)
handles = [mpatches.Patch(color=LOC_C[l], label=l.title()) for l in LOC_ORDER]
ax_pca.legend(handles=handles, fontsize=10, frameon=True, loc='upper right')

# ── Panel B: PERMANOVA R² bar chart for multiple variables ────────────────────
ax_perm = fig.add_subplot(gs[0, 2])

X_hel = hellinger(DATA.values)
D_all = bray_curtis(X_hel)

perm_vars = {
    'Anatomical site':    meta['Location'].values,
    'Ig fraction':        meta['Antibody'].values,
    'Clinical condition': meta['Condition'].values,
    'Storage method':     meta['Storage'].values,
}
perm_res = {}
for label, grp in perm_vars.items():
    perm_res[label] = permanova(D_all, grp, n_perm=499)

labels_p = list(perm_res.keys())
r2_vals  = [perm_res[l]['R2'] for l in labels_p]
p_vals   = [perm_res[l]['p']  for l in labels_p]
bar_cols  = ['#2166AC', '#E66101', '#D7191C', '#1565C0']

bars = ax_perm.barh(range(len(labels_p)), r2_vals, color=bar_cols, alpha=0.80,
                    edgecolor='white', height=0.65)
for i, (r2, p) in enumerate(zip(r2_vals, p_vals)):
    ax_perm.text(r2 + 0.002, i, f'R²={r2:.3f} {pct_label(p)}',
                 va='center', fontsize=10, fontweight='bold')

ax_perm.set_yticks(range(len(labels_p)))
ax_perm.set_yticklabels(labels_p, fontsize=11)
ax_perm.set_xlabel('PERMANOVA R² (variance explained)', fontsize=10)
ax_perm.set_title('B   PERMANOVA effect sizes\n(Bray-Curtis, n=499 permutations)',
                  fontsize=11, fontweight='bold', loc='left', pad=6)
ax_perm.spines[['top', 'right']].set_visible(False)
ax_perm.xaxis.grid(True, alpha=0.22, linestyle='--')
ax_perm.set_axisbelow(True)
ax_perm.set_xlim(0, max(r2_vals) * 1.4)

# ── Panel C: Spatial Ig gradient for key taxa ─────────────────────────────────
ax_grad = fig.add_subplot(gs[1, :])

KEY_TAXA_GRAD = [
    'Lactobacillus iners', 'Lactobacillus crispatus', 'Lactobacillus gasseri',
    'Gardnerella vaginalis', 'Sneathia vaginalis', 'Fannyhessea vaginae',
    'Parvimonas micra', 'Prevotella buccalis'
]
KEY_TAXA_GRAD = [t for t in KEY_TAXA_GRAD if t in taxa]
spp_palette = plt.cm.tab10(np.linspace(0, 1, max(1, len(KEY_TAXA_GRAD))))

x_pos    = 0
xtick_pos, xtick_lab = [], []
n_ab     = len(['A', 'M', 'G'])
spp_centers, spp_labels, spp_colors = [], [], []
xtick_colors = []

for ti, taxon in enumerate(KEY_TAXA_GRAD):
    group_positions = []
    spp_color = spp_palette[ti % len(spp_palette)]
    for ab, ab_label in [('A', 'IgA'), ('M', 'IgM'), ('G', 'IgG')]:
        vals_per_loc = []
        for loc in LOC_ORDER:
            sids = meta[(meta['Location'] == loc) & (meta['Antibody'] == ab)].index
            sids = [s for s in sids if s in DATA.index]
            mu   = DATA.loc[sids, taxon].mean() * 100 if sids else 0
            vals_per_loc.append(mu)

        col = AB_C[ab]
        ax_grad.bar(x_pos, vals_per_loc[0], width=0.25, color=LOC_C['VAGINA'],
                    alpha=0.85, edgecolor='white', linewidth=0.3)
        ax_grad.bar(x_pos + 0.25, vals_per_loc[1], width=0.25, color=LOC_C['CERVIX'],
                    alpha=0.85, edgecolor='white', linewidth=0.3)
        ax_grad.bar(x_pos + 0.50, vals_per_loc[2], width=0.25,
                    color=LOC_C['ENDOMETRIUM'],
                    alpha=0.85, edgecolor='white', linewidth=0.3)
        tick_center = x_pos + 0.25
        group_positions.append(tick_center)
        xtick_pos.append(tick_center)
        xtick_lab.append(ab_label)
        xtick_colors.append(spp_color)
        x_pos += 0.9

    t_short = taxon.replace('Lactobacillus ', 'L. ').replace('Gardnerella ', 'G. ')
    spp_centers.append(float(np.mean(group_positions)))
    spp_labels.append(t_short)
    spp_colors.append(spp_color)
    x_pos += 0.5

ax_grad.set_xlim(-0.5, x_pos)
ax_grad.set_ylabel('Mean relative abundance (%)', fontsize=11)
ax_grad.set_title(
    'C   Spatial Ig gradient — key taxon abundance across Ig fractions and anatomical sites (bar groups: Vagina / Cervix / Endometrium)',
    fontsize=11, fontweight='bold', loc='left', pad=10)
ax_grad.spines[['top', 'right']].set_visible(False)
ax_grad.yaxis.grid(True, alpha=0.22, linestyle='--')
ax_grad.set_axisbelow(True)
ax_grad.set_xticks(xtick_pos)
ax_grad.set_xticklabels(xtick_lab, fontsize=8.5, fontweight='bold')
ax_grad.tick_params(axis='x', length=0, pad=2)
for tick, tick_color in zip(ax_grad.get_xticklabels(), xtick_colors):
    tick.set_color(tick_color)

if spp_centers:
    secax = ax_grad.secondary_xaxis('bottom')
    secax.set_xticks(spp_centers)
    secax.set_xticklabels(spp_labels, fontsize=8, fontstyle='italic', fontweight='bold')
    secax.tick_params(axis='x', length=0, pad=18)
    secax.spines['bottom'].set_visible(False)
    for lbl, s_color in zip(secax.get_xticklabels(), spp_colors):
        lbl.set_color(s_color)

# Colour legend
loc_handles = [mpatches.Patch(color=LOC_C[l], label=l.title()) for l in LOC_ORDER]
ab_handles  = [mpatches.Patch(color=AB_C[ab], alpha=0.7, label=AB_LABELS[ab])
               for ab in ['A', 'M', 'G']]
ax_grad.legend(handles=loc_handles, fontsize=9, frameon=True,
               loc='upper right', ncol=6)

plt.tight_layout()
fig.savefig(os.path.join(FIG_DIR, 'Figure9.png'),
            dpi=600, bbox_inches='tight', facecolor='white')
plt.close()
print("Saved Figure9.png")
