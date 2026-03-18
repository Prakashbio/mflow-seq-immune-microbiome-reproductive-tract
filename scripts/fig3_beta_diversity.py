import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Ellipse
import warnings; warnings.filterwarnings('ignore')

from shared_utils import (load_data, hellinger, bray_curtis, pcoa, nmds,
                           permanova, pct_label,
                           LOC_C, AB_C, COND_C, STOR_C,
                           FIG_DIR)

otu_rel, meta, tax = load_data()

# ── Hellinger transform then Bray-Curtis ──────────────────────────────────────
X        = otu_rel.T.values        # (n_samples × n_taxa)
samples  = otu_rel.columns.tolist()
X_hel    = hellinger(X)
D        = bray_curtis(X_hel)

# ── Ordinations ───────────────────────────────────────────────────────────────
coords_pcoa, var_pcoa = pcoa(D, n_axes=4)

coords_nmds = nmds(D, n_axes=2, n_init=5, max_iter=200)

# ── PERMANOVA for each grouping variable ─────────────────────────────────────
perm_results = {}
for col, label in [('Location','Location'), ('Antibody','Antibody fraction'),
                    ('Condition','Clinical condition'), ('Storage','Storage method')]:
    grp = meta.loc[samples, col].values
    perm_results[label] = permanova(D, grp, n_perm=499)

# ── Confidence ellipse helper ─────────────────────────────────────────────────
def plot_ellipse(ax, pts, color, alpha=0.18, n_std=1.96):
    if len(pts) < 3:
        return
    mean  = pts.mean(axis=0)
    cov   = np.cov(pts.T)
    if cov.ndim < 2 or np.any(np.isnan(cov)):
        return
    vals, vecs = np.linalg.eigh(cov)
    order  = vals.argsort()[::-1]
    vals   = vals[order]; vecs = vecs[:, order]
    angle  = np.degrees(np.arctan2(*vecs[:, 0][::-1]))
    w, h   = 2 * n_std * np.sqrt(np.abs(vals))
    ell    = Ellipse(xy=mean, width=w, height=h, angle=angle,
                     facecolor=color, alpha=alpha, edgecolor=color, linewidth=1.5)
    ax.add_patch(ell)


# ── Build colour mappings ─────────────────────────────────────────────────────
COLOR_MAPS = {
    'Location':          (meta.loc[samples, 'Location'].values,  LOC_C,  'Location'),
    'Antibody fraction': (meta.loc[samples, 'Antibody'].values,  AB_C,   'Antibody fraction'),
    'Clinical condition':(meta.loc[samples, 'Condition'].values, COND_C, 'Condition'),
    'Storage method':    (meta.loc[samples, 'Storage'].values,   STOR_C, 'Storage'),
}


def draw_ordination(ax, coords, labels, color_map, title, xlabel, ylabel):
    unique_labs = np.unique(labels)
    for lab in unique_labs:
        idx  = np.where(labels == lab)[0]
        pts  = coords[idx, :2]
        col  = color_map.get(lab, '#888888')
        ax.scatter(pts[:, 0], pts[:, 1], c=col, s=55, alpha=0.82,
                   edgecolors='white', lw=0.4, label=lab.title(), zorder=5)
        plot_ellipse(ax, pts, col)
    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_ylabel(ylabel, fontsize=9)
    ax.set_title(title, fontsize=11, fontweight='bold', pad=5)
    ax.spines[['top', 'right']].set_visible(False)
    ax.axhline(0, color='#ccc', lw=0.8, linestyle='--')
    ax.axvline(0, color='#ccc', lw=0.8, linestyle='--')


# ─────────────────────────────────────────────────────────────────────────────
# MAIN FIGURE
# ─────────────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 4, figsize=(24, 12))
fig.patch.set_facecolor('white')

col_keys = ['Location', 'Antibody fraction', 'Clinical condition', 'Storage method']
pcoord_xy = (f"PCoA1 ({var_pcoa[0]*100:.1f}%)", f"PCoA2 ({var_pcoa[1]*100:.1f}%)")

for ci, key in enumerate(col_keys):
    labels, cmap, cmap_title = COLOR_MAPS[key]

    # PCoA row
    draw_ordination(axes[0, ci], coords_pcoa, labels, cmap,
                    f'PCoA — {key}', pcoord_xy[0], pcoord_xy[1] if ci == 0 else '')
    pr = perm_results[key]
    axes[0, ci].text(0.97, 0.03,
                     f"PERMANOVA\nR²={pr['R2']:.3f}, p={pr['p']:.3f} {pct_label(pr['p'])}",
                     transform=axes[0, ci].transAxes, ha='right', va='bottom',
                     fontsize=8, color='#333', style='italic',
                     bbox=dict(boxstyle='round,pad=0.3', fc='#f8f8f8', ec='#ccc', alpha=0.85))

    # NMDS row
    draw_ordination(axes[1, ci], coords_nmds, labels, cmap,
                    f'NMDS — {key}', 'NMDS1', 'NMDS2' if ci == 0 else '')

    # Add legend
    unique_labs = np.unique(labels)
    handles = [mpatches.Patch(color=cmap.get(l, '#888'), alpha=0.8,
                              label=l.title()) for l in unique_labs]
    axes[1, ci].legend(handles=handles, fontsize=8, frameon=True,
                       loc='lower right', markerscale=1.2)

plt.tight_layout()
fig.savefig(os.path.join(FIG_DIR, 'Figure3.png'),
            dpi=600, bbox_inches='tight', facecolor='white')
plt.close()
print("Saved Figure3.png")
