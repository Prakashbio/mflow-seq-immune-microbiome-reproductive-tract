import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from scipy.stats import spearmanr
import warnings; warnings.filterwarnings('ignore')

from shared_utils import (load_data, LOC_C, AB_C, LOC_ORDER, AB_ORDER, AB_LABELS, FIG_DIR)

otu_rel, meta, tax = load_data()
DATA  = otu_rel.T
taxa  = otu_rel.index.tolist()


def spread_positions(values, min_gap=0.38):
    if not values:
        return []
    order = np.argsort(values)
    adjusted = np.array(values, dtype=float)[order]
    for i in range(1, len(adjusted)):
        adjusted[i] = max(adjusted[i], adjusted[i - 1] + min_gap)
    for i in range(len(adjusted) - 2, -1, -1):
        adjusted[i] = min(adjusted[i], adjusted[i + 1] - min_gap)
    result = np.empty_like(adjusted)
    result[np.argsort(order)] = adjusted
    return result.tolist()

# ─────────────────────────────────────────────────────────────────────────────
# MAIN FIGURE
# ─────────────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(22, 16))
fig.patch.set_facecolor('white')

gs = fig.add_gridspec(2, 2, hspace=0.15, wspace=0.38)

# ── Panel A: Filtering summary ────────────────────────────────────────────────
ax_filt = fig.add_subplot(gs[0, 0])
filt = pd.read_excel(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                   '..', 'data', 'filtering_statistics.xlsx'),
                     header=0)
group_parts = filt['Group'].str.rsplit('_', n=1, expand=True)
filt['Location'] = group_parts[0]
filt['Antibody'] = group_parts[1]
filt = filt[filt['Antibody'].isin(['A', 'G', 'M'])].copy()

thresh_cols = [c for c in filt.columns if 'Bacteria_Count' in c]
thresh_labs = ['>1% ≥2 samples', '>1% ≥3 samples', '>1% ≥4 samples', '>10% ≥2 samples']
ab_order_panel = ['A', 'G', 'M']
ab_label_panel = {'A': 'IgA', 'G': 'IgG', 'M': 'IgM'}

x_pos = 0.0
x, xtick_labels, xtick_abs, xtick_locs = [], [], [], []
organ_centers, organ_labels, organ_keys = [], [], []
plot_rows = []
for loc in LOC_ORDER:
    group_x = []
    for ab in ab_order_panel:
        row = filt[(filt['Location'] == loc) & (filt['Antibody'] == ab)]
        if row.empty:
            continue
        plot_rows.append(row.iloc[0])
        x.append(x_pos)
        xtick_labels.append(ab_label_panel[ab])
        xtick_abs.append(ab)
        xtick_locs.append(loc)
        group_x.append(x_pos)
        x_pos += 1.0
    if group_x:
        organ_centers.append(float(np.mean(group_x)))
        organ_labels.append(loc.title())
        organ_keys.append(loc)
        x_pos += 0.55

plot_df = pd.DataFrame(plot_rows).reset_index(drop=True)
n_samples = plot_df['Number of Samples in Group'].values

bw  = 0.18
colors_filt = ['#4393C3', '#2166AC', '#08306B', '#E66101']

for i, (col, lab, col_c) in enumerate(zip(thresh_cols, thresh_labs, colors_filt)):
    ax_filt.bar(np.array(x) + i * bw, plot_df[col].values, width=bw, label=lab,
                color=col_c, alpha=0.82, edgecolor='white')

tick_pos = np.array(x) + bw * 1.5
ax_filt.set_xticks(tick_pos)
ax_filt.set_xticklabels(xtick_labels, fontsize=8.5, fontweight='bold')
for label, loc in zip(ax_filt.get_xticklabels(), xtick_locs):
    label.set_color(LOC_C[loc])
if organ_centers:
    secax = ax_filt.secondary_xaxis('bottom')
    secax.set_xticks(organ_centers)
    secax.set_xticklabels(organ_labels, fontsize=8.5, fontweight='bold')
    secax.tick_params(axis='x', length=0, pad=18)
    secax.spines['bottom'].set_visible(False)
    for label, loc in zip(secax.get_xticklabels(), organ_keys):
        label.set_color(LOC_C[loc])
ax_filt.set_ylabel('Number of taxa retained', fontsize=10)
ax_filt.set_title('A   Filtering threshold comparison', fontsize=11,
                  fontweight='bold', loc='left', pad=6)
ax_filt.spines[['top', 'right']].set_visible(False)
ax_filt.yaxis.grid(True, alpha=0.22, linestyle='--')
ax_filt.set_axisbelow(True)
if len(tick_pos):
    ax_filt.set_xlim(tick_pos.min() - 0.5, tick_pos.max() + 0.5)

# Add sample size as secondary info
ax2 = ax_filt.twinx()
sample_line = ax2.plot(tick_pos, n_samples, 'k--o', ms=5, lw=1.5, alpha=0.55, label='N samples')[0]
ax2.set_ylabel('N samples in group', fontsize=9, color='#555')
ax2.tick_params(axis='y', labelcolor='#555')
ax2.spines[['top']].set_visible(False)

# ── Panel B: Prevalence heatmap ───────────────────────────────────────────────
ax_prev = fig.add_subplot(gs[0, 1])

prev_rows = []
prev_abs = []
prev_loc_keys = []
prev_org_centers, prev_org_labels, prev_org_keys = [], [], []
prev_ab_order = AB_ORDER
prev_ab_labels = {ab: AB_LABELS[ab] for ab in AB_ORDER}
prev_col_idx = 0
for loc in LOC_ORDER:
    group_cols = []
    for ab in prev_ab_order:
        sids = meta[(meta['Location'] == loc) & (meta['Antibody'] == ab)].index
        sids = [s for s in sids if s in DATA.index]
        if not sids: continue
        sub = DATA.loc[sids]
        # prevalence = fraction of samples where abundance > 1%
        prev = (sub > 0.01).mean(axis=0).values
        prev_rows.append(prev)
        prev_abs.append(ab)
        prev_loc_keys.append(loc)
        group_cols.append(prev_col_idx)
        prev_col_idx += 1
    if group_cols:
        prev_org_centers.append(float(np.mean(group_cols)))
        prev_org_labels.append(loc.title())
        prev_org_keys.append(loc)

prev_mat = np.array(prev_rows).T  # taxa × groups

# Sort by total prevalence
sort_idx  = prev_mat.sum(axis=1).argsort()[::-1]
prev_mat  = prev_mat[sort_idx]
taxa_prev = [taxa[i] for i in sort_idx]

im = ax_prev.imshow(prev_mat, aspect='auto', cmap='Blues', vmin=0, vmax=1,
                     interpolation='nearest')
ax_prev.set_xticks(range(len(prev_abs)))
ax_prev.set_xticklabels([prev_ab_labels[ab] for ab in prev_abs], fontsize=8.5, fontweight='bold')
for label, loc in zip(ax_prev.get_xticklabels(), prev_loc_keys):
    label.set_color(LOC_C[loc])
if prev_org_centers:
    prev_secax = ax_prev.secondary_xaxis('bottom')
    prev_secax.set_xticks(prev_org_centers)
    prev_secax.set_xticklabels(prev_org_labels, fontsize=8.5, fontweight='bold')
    prev_secax.tick_params(axis='x', length=0, pad=18)
    prev_secax.spines['bottom'].set_visible(False)
    for label, loc in zip(prev_secax.get_xticklabels(), prev_org_keys):
        label.set_color(LOC_C[loc])
ax_prev.set_yticks(range(len(taxa_prev)))
ax_prev.set_yticklabels(taxa_prev, fontsize=7.5, style='italic')
ax_prev.set_title('B   Taxon prevalence (fraction of samples >1%)',
                  fontsize=11, fontweight='bold', loc='left', pad=6)

# Organ group dividers
n_ab = len(prev_ab_order)
for xi in [n_ab - 0.5, 2 * n_ab - 0.5]:
    ax_prev.axvline(xi, color='white', lw=3)

plt.colorbar(im, ax=ax_prev, label='Prevalence (prop. samples >1%)',
             fraction=0.03, pad=0.01)

# ── Panel C: Co-occurrence network (Pre-sorted all organs) ────────────────────
ax_net = fig.add_subplot(gs[1, :])

pre_sids = meta[meta['Antibody'] == 'PRE'].index
pre_sids = [s for s in pre_sids if s in DATA.index]
pre_data = DATA.loc[pre_sids]

# Spearman correlation between taxa
corr_mat = np.zeros((len(taxa), len(taxa)))
pval_mat = np.ones((len(taxa), len(taxa)))
for i in range(len(taxa)):
    for j in range(i, len(taxa)):
        x_v = pre_data.iloc[:, i].values
        y_v = pre_data.iloc[:, j].values
        if np.std(x_v) > 0 and np.std(y_v) > 0:
            r, p = spearmanr(x_v, y_v)
            corr_mat[i, j] = corr_mat[j, i] = r
            pval_mat[i, j] = pval_mat[j, i] = p
        if i == j:
            corr_mat[i, j] = 1.0

# Network: only significant correlations (p < 0.05, |r| > 0.4)
THRESH_R = 0.4
THRESH_P = 0.05

# Node positions using simple circular layout
n_taxa = len(taxa)
angles = np.linspace(0, 2 * np.pi, n_taxa, endpoint=False)
radius = 5.0
node_x = radius * np.cos(angles)
node_y = radius * np.sin(angles)

# Node size = mean abundance
node_sizes = pre_data.mean(axis=0).values * 2000 + 40

# Colour nodes by phylum — align taxonomy to current taxa (post-decontamination)
phyla = [tax.loc[t, 'Phylum'] if t in tax.index else 'Unknown' for t in taxa]
phylum_colors_net = {
    'Bacillota':      '#4393C3',
    'Actinomycetota': '#D6604D',
    'Pseudomonadota': '#74C476',
    'Bacteroidota':   '#FD8D3C',
    'Fusobacteriota': '#9E9AC8',
    'Mycoplasmatota': '#FDAE6B',
}
node_colors_net = [phylum_colors_net.get(p, '#BBBBBB') for p in phyla]

# Draw edges
edge_count = 0
for i in range(n_taxa):
    for j in range(i + 1, n_taxa):
        r = corr_mat[i, j]; p = pval_mat[i, j]
        if abs(r) >= THRESH_R and p < THRESH_P:
            edge_col = '#E66101' if r > 0 else '#5E3C99'
            lw = min(abs(r) * 4, 4)
            ax_net.plot([node_x[i], node_x[j]], [node_y[i], node_y[j]],
                        color=edge_col, lw=lw, alpha=0.45, zorder=1)
            edge_count += 1

# Draw nodes
ax_net.scatter(node_x, node_y, s=node_sizes, c=node_colors_net, alpha=0.90,
               edgecolors='white', lw=1.0, zorder=5)
# Labels
label_specs = []
for tx, x_n, y_n in zip(taxa, node_x, node_y):
    side = 'right' if x_n >= 0 else 'left'
    label_specs.append({
        'taxon': tx,
        'x_node': x_n,
        'y_node': y_n,
        'side': side,
        'x_text': x_n * 1.18,
        'y_text': y_n * 1.18,
    })

for side in ['left', 'right']:
    side_specs = [spec for spec in label_specs if spec['side'] == side]
    if not side_specs:
        continue
    adjusted_y = spread_positions([spec['y_text'] for spec in side_specs], min_gap=0.20)
    for spec, new_y in zip(side_specs, adjusted_y):
        spec['y_text'] = new_y

label_x = []
label_y = []
for spec in label_specs:
    label_x.append(spec['x_text'])
    label_y.append(spec['y_text'])
    ax_net.annotate(spec['taxon'], xy=(spec['x_node'], spec['y_node']),
                    xytext=(spec['x_text'], spec['y_text']),
                    fontsize=7.5, ha='left' if spec['side'] == 'right' else 'right',
                    va='center', style='italic', color='#222')

x_all = np.concatenate([node_x, np.array(label_x)])
y_all = np.concatenate([node_y, np.array(label_y)])
ax_net.set_xlim(x_all.min() - 1.0, x_all.max() + 1.0)
ax_net.set_ylim(y_all.min() - 1.0, y_all.max() + 0.35)
ax_net.set_aspect('equal')
ax_net.axis('off')
ax_net.set_title(
    f'C   Microbial co-occurrence network (Pre-sorted; |Spearman r| ≥ {THRESH_R}, p < {THRESH_P};  '
    f'{edge_count} edges shown)',
    fontsize=11, fontweight='bold', loc='left', pad=0)

# Network legend
edge_handles = [Line2D([0], [0], color='#E66101', lw=2.5, label='Positive correlation'),
                Line2D([0], [0], color='#5E3C99', lw=2.5, label='Negative correlation')]
phylum_handles = [mpatches.Patch(color=c, alpha=0.8, label=p)
                  for p, c in phylum_colors_net.items()]
ax_net.legend(handles=edge_handles + phylum_handles, fontsize=8.5,
              frameon=True, loc='lower right', bbox_to_anchor=(1.0, -0.03), ncol=4)

plt.tight_layout()

# Dedicated legend axis below Panel A.
filt_handles = [mpatches.Patch(facecolor=c, alpha=0.82, edgecolor='white', label=l)
                for c, l in zip(colors_filt, thresh_labs)] + [sample_line]
filt_pos = ax_filt.get_position()
prev_pos = ax_prev.get_position()
net_pos = ax_net.get_position()
leg_h = 0.03
leg_y = max(net_pos.y1 + 0.002, min(filt_pos.y0, prev_pos.y0) - 0.045)
legend_ax = fig.add_axes([filt_pos.x0, leg_y, filt_pos.width, leg_h])
legend_ax.axis('off')
legend_ax.legend(handles=filt_handles, loc='center', ncol=5, frameon=False,
                 fontsize=8, columnspacing=1.4, handlelength=1.8)

fig.savefig(os.path.join(FIG_DIR, 'SupplementaryFigure14.png'),
            dpi=600, bbox_inches='tight', facecolor='white')
plt.close()
print("Saved SupplementaryFigure14.png")
