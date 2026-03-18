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

from shared_utils import (load_data, shannon, hellinger, bray_curtis, pct_label,
                           LOC_C, COND_C, LOC_ORDER, AB_ORDER, AB_LABELS, FIG_DIR)

# Load raw counts too
otu_raw = pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    '..', 'data', 'otu_table.csv'), index_col=0)
otu_rel, meta, tax = load_data()
DATA  = otu_rel.T
taxa  = otu_rel.index.tolist()
meta['Patient'] = meta.index.str.split('_').str[0]

# ─────────────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(22, 18))
fig.patch.set_facecolor('white')
gs = fig.add_gridspec(2, 3, hspace=0.15, wspace=0.15)

# ── Panel A: Per-patient spaghetti plot (pre-sorted, Shannon by organ) ────────
ax_a = fig.add_subplot(gs[0, 0:2])

patients_all3 = []
for pat in meta['Patient'].unique():
    locs = meta[(meta['Patient'] == pat) & (meta['Antibody'] == 'PRE')]['Location'].unique()
    if set(LOC_ORDER).issubset(set(locs)):
        patients_all3.append(pat)

print(f"Patients with all 3 organs (PRE fraction): {len(patients_all3)}")

cond_per_pat = {}
for pat in patients_all3:
    sid = meta[(meta['Patient'] == pat) & (meta['Antibody'] == 'PRE') &
               (meta['Location'] == 'VAGINA')].index
    if len(sid):
        cond_per_pat[pat] = meta.loc[sid[0], 'Condition']

x_pos = {'VAGINA': 0, 'CERVIX': 1, 'ENDOMETRIUM': 2}

for pat in patients_all3:
    cond = cond_per_pat.get(pat, 'UNKNOWN')
    col  = COND_C.get(cond, '#888888')
    pts  = {}
    for loc in LOC_ORDER:
        sid = meta[(meta['Patient'] == pat) & (meta['Antibody'] == 'PRE') &
                   (meta['Location'] == loc)].index
        sid = [s for s in sid if s in DATA.index]
        if sid:
            pts[loc] = shannon(DATA.loc[sid[0]].values)
    if len(pts) == 3:
        xs = [x_pos[l] for l in LOC_ORDER if l in pts]
        ys = [pts[l] for l in LOC_ORDER if l in pts]
        ax_a.plot(xs, ys, '-o', color=col, alpha=0.55, lw=1.8, ms=7,
                  markeredgecolor='white', markeredgewidth=0.5, zorder=4)
        ax_a.text(2.05, ys[-1], pat, fontsize=7.5, va='center', color=col, alpha=0.7)

# Group means
for loc in LOC_ORDER:
    for cond in ['HEALTHY', 'DYSBIOTIC']:
        vals = []
        for pat in patients_all3:
            if cond_per_pat.get(pat) != cond: continue
            sid = meta[(meta['Patient'] == pat) & (meta['Antibody'] == 'PRE') &
                       (meta['Location'] == loc)].index
            sid = [s for s in sid if s in DATA.index]
            if sid: vals.append(shannon(DATA.loc[sid[0]].values))
        if vals:
            xp = x_pos[loc] + (0.08 if cond == 'HEALTHY' else -0.08)
            ax_a.scatter(xp, np.mean(vals), s=180, color=COND_C[cond],
                         zorder=10, edgecolors='white', lw=1.5,
                         marker='D', label=f'{cond.title()} mean' if loc == 'VAGINA' else '_')

ax_a.set_xticks([0, 1, 2])
ax_a.set_xticklabels(['Vagina', 'Cervix', 'Endometrium'], fontsize=12, fontweight='bold')
for tick, loc in zip(ax_a.get_xticklabels(), LOC_ORDER):
    tick.set_color(LOC_C[loc])
ax_a.set_ylabel("Shannon diversity H'", fontsize=11)
ax_a.set_title(f'A   Per-patient alpha diversity trajectories (Pre-sorted, n={len(patients_all3)} with all 3 sites)',
               fontsize=11, fontweight='bold', loc='left', pad=6)
ax_a.spines[['top', 'right']].set_visible(False)
ax_a.yaxis.grid(True, alpha=0.22, linestyle='--')
ax_a.set_axisbelow(True)
handles_a = [mpatches.Patch(color=COND_C['HEALTHY'],  alpha=0.6, label='Healthy'),
             mpatches.Patch(color=COND_C['DYSBIOTIC'], alpha=0.6, label='Dysbiotic'),
             Line2D([0], [0], marker='D', color='w', markerfacecolor='#555',
                    markersize=10, label='Group mean')]
ax_a.legend(handles=handles_a, fontsize=10, frameon=True, loc='lower right')

# ── Panel B: Read depth distribution ────────────────────────────────────────
ax_b = fig.add_subplot(gs[0, 2])

read_depths = otu_raw.sum(axis=0)
rdf = pd.DataFrame({'ReadDepth': read_depths, 'SampleID': read_depths.index})
rdf['Location'] = rdf['SampleID'].apply(lambda x: meta.loc[x, 'Location'] if x in meta.index else 'UNKNOWN')
rdf['Antibody'] = rdf['SampleID'].apply(lambda x: meta.loc[x, 'Antibody'] if x in meta.index else 'UNKNOWN')

for ai, ab in enumerate(AB_ORDER):
    sub_rd = rdf[rdf['Antibody'] == ab]['ReadDepth'].values
    if len(sub_rd) == 0: continue
    col = {'PRE': '#636363', 'A': '#E66101', 'M': '#5E3C99', 'G': '#1B9E77'}[ab]
    bp = ax_b.boxplot(sub_rd, positions=[ai], widths=0.6, patch_artist=True,
                      showfliers=True,
                      medianprops={'color': 'white', 'lw': 2.5},
                      whiskerprops={'color': col, 'lw': 1.5},
                      capprops={'color': col, 'lw': 1.5},
                      flierprops={'marker': 'o', 'ms': 4, 'alpha': 0.5,
                                   'markerfacecolor': col},
                      boxprops={'facecolor': col, 'alpha': 0.55,
                                 'edgecolor': col, 'lw': 1.2})

ax_b.axhline(10000, color='#E31A1C', lw=1.5, linestyle='--', alpha=0.7,
             label='10,000 reads threshold')
ax_b.set_yscale('log')
ax_b.set_xticks(range(len(AB_ORDER)))
ax_b.set_xticklabels([AB_LABELS[ab] for ab in AB_ORDER], fontsize=11, fontweight='bold')
for tick, ab in zip(ax_b.get_xticklabels(), AB_ORDER):
    tick.set_color({'PRE': '#636363', 'A': '#E66101', 'M': '#5E3C99', 'G': '#1B9E77'}[ab])
ax_b.set_ylabel('Sequencing read depth (log scale)', fontsize=10)
ax_b.set_title('B   Read depth per Ig fraction\n(quality control)', fontsize=11,
               fontweight='bold', loc='left', pad=6)
ax_b.spines[['top', 'right']].set_visible(False)
ax_b.yaxis.grid(True, alpha=0.22, linestyle='--')
ax_b.set_axisbelow(True)
ax_b.legend(fontsize=9, frameon=True, loc='upper right')
ax_b.text(0.03, 0.03,
          f"n samples with <5k reads: {(read_depths < 5000).sum()}\n"
          f"min: {read_depths.min():.0f}  max: {read_depths.max():.0f}",
          transform=ax_b.transAxes, fontsize=8.5, va='bottom',
          bbox=dict(boxstyle='round,pad=0.3', fc='#f8f8f8', ec='#ccc', alpha=0.85))

# ── Panel C: Within-patient vs between-patient Bray-Curtis ───────────────────
ax_c = fig.add_subplot(gs[1, 0:2])

# For each pair of anatomical sites, compute within-patient BC distances
# (same patient, different organs) vs between-patient (different patients)
site_pairs = [('VAGINA', 'CERVIX'), ('VAGINA', 'ENDOMETRIUM'), ('CERVIX', 'ENDOMETRIUM')]
pair_labels = ['Vagina\nvs Cervix', 'Vagina\nvs\nEndometrium', 'Cervix\nvs\nEndometrium']
within_all  = {p: [] for p in pair_labels}
between_all = {p: [] for p in pair_labels}

for (loc1, loc2), plabel in zip(site_pairs, pair_labels):
    pats_with_both = []
    for pat in meta['Patient'].unique():
        s1 = meta[(meta['Patient']==pat)&(meta['Location']==loc1)&(meta['Antibody']=='PRE')].index
        s2 = meta[(meta['Patient']==pat)&(meta['Location']==loc2)&(meta['Antibody']=='PRE')].index
        s1 = [s for s in s1 if s in DATA.index]
        s2 = [s for s in s2 if s in DATA.index]
        if s1 and s2:
            pats_with_both.append((pat, s1[0], s2[0]))

    for i, (pat_i, s1_i, s2_i) in enumerate(pats_with_both):
        v1 = hellinger(DATA.loc[s1_i].values.reshape(1, -1))
        v2 = hellinger(DATA.loc[s2_i].values.reshape(1, -1))
        from scipy.spatial.distance import braycurtis
        within_all[plabel].append(braycurtis(v1.flatten(), v2.flatten()))
        # Between: same site1 but different patient's site2
        for j, (pat_j, s1_j, s2_j) in enumerate(pats_with_both):
            if i >= j: continue
            v1j = hellinger(DATA.loc[s1_j].values.reshape(1, -1))
            v2j = hellinger(DATA.loc[s2_j].values.reshape(1, -1))
            between_all[plabel].append(braycurtis(v1.flatten(), v2j.flatten()))

x_pos_c = 0
xtick_c, xtick_c_lab = [], []
for pi, plabel in enumerate(pair_labels):
    w = within_all[plabel]
    b = between_all[plabel]
    for xi, (vals, label, col, hatch) in enumerate([
        (w, 'Within-patient', '#2166AC', ''),
        (b, 'Between-patient', '#D7191C', '///')
    ]):
        if not vals: continue
        bp = ax_c.boxplot(vals, positions=[x_pos_c], widths=0.7, patch_artist=True,
                          showfliers=True,
                          medianprops={'color': 'white', 'lw': 2.5},
                          whiskerprops={'color': col, 'lw': 1.5},
                          capprops={'color': col, 'lw': 1.5},
                          flierprops={'marker': 'o', 'ms': 3, 'alpha': 0.4,
                                       'markerfacecolor': col},
                          boxprops={'facecolor': col, 'alpha': 0.50,
                                     'edgecolor': col, 'lw': 1.2,
                                     'hatch': hatch})
        xtick_c.append(x_pos_c)
        xtick_c_lab.append(f'{label}\n({plabel.replace(chr(10)," ")})')
        x_pos_c += 1

    # MWU
    if within_all[plabel] and between_all[plabel]:
        try:
            _, p = mannwhitneyu(within_all[plabel], between_all[plabel], alternative='two-sided')
            xm = x_pos_c - 1.5
            ymax = max(max(within_all[plabel]), max(between_all[plabel]))
            ax_c.annotate('', xy=(x_pos_c - 1, ymax + 0.03),
                          xytext=(x_pos_c - 2, ymax + 0.03),
                          arrowprops=dict(arrowstyle='-', color='#555', lw=1.5))
            ax_c.text(x_pos_c - 1.5, ymax + 0.05, pct_label(p),
                      ha='center', fontsize=11, fontweight='bold')
        except: pass
    x_pos_c += 0.6

ax_c.set_xticks(xtick_c)
ax_c.set_xticklabels(xtick_c_lab, fontsize=8.5)
ax_c.set_ylabel('Bray-Curtis dissimilarity', fontsize=10)
ax_c.set_title('C   Within- vs between-patient community similarity across anatomical sites (Pre-sorted fraction, paired analysis)',
               fontsize=11, fontweight='bold', loc='left', pad=15)
ax_c.set_ylim(bottom=0)
ax_c.spines[['top', 'right']].set_visible(False)
ax_c.yaxis.grid(True, alpha=0.22, linestyle='--')
ax_c.set_axisbelow(True)
handles_c = [mpatches.Patch(color='#2166AC', alpha=0.6, label='Within-patient (same individual, different sites)'),
             mpatches.Patch(color='#D7191C', alpha=0.6, label='Between-patient (different individuals)')]
ax_c.legend(handles=handles_c, fontsize=9, frameon=True, loc='lower right')

# ── Panel D: Read depth vs Shannon scatter ───────────────────────────────────
ax_d = fig.add_subplot(gs[1, 2])

for loc in LOC_ORDER:
    sids = [s for s in meta[meta['Location'] == loc].index if s in DATA.index]
    rd   = [otu_raw[s].sum() for s in sids]
    sh   = [shannon(DATA.loc[s].values) for s in sids]
    ax_d.scatter(rd, sh, c=LOC_C[loc], s=45, alpha=0.70,
                 edgecolors='white', lw=0.4, label=loc.title(), zorder=5)

# Spearman correlation
all_sids = [s for s in meta.index if s in DATA.index]
all_rd   = np.array([otu_raw[s].sum() for s in all_sids])
all_sh   = np.array([shannon(DATA.loc[s].values) for s in all_sids])
from scipy.stats import spearmanr
r_sp, p_sp = spearmanr(np.log10(all_rd + 1), all_sh)
ax_d.text(0.05, 0.95,
          f'Spearman r = {r_sp:.3f}\np = {p_sp:.4f}',
          transform=ax_d.transAxes, fontsize=10,
          bbox=dict(boxstyle='round,pad=0.3', fc='#f8f8f8', ec='#ccc', alpha=0.85))
ax_d.axvline(5000, color='#E31A1C', lw=1.5, linestyle='--', alpha=0.7,
             label='5,000 reads')
ax_d.set_xscale('log')
ax_d.set_xlabel('Sequencing read depth (log scale)', fontsize=10)
ax_d.set_ylabel("Shannon H'", fontsize=10)
ax_d.set_title('D   Read depth vs Shannon diversity\n(correlation check for depth bias)',
               fontsize=11, fontweight='bold', loc='left', pad=6)
ax_d.spines[['top', 'right']].set_visible(False)
ax_d.grid(True, alpha=0.18, linestyle='--')
handles_d = [mpatches.Patch(color=LOC_C[l], label=l.title()) for l in LOC_ORDER]
ax_d.legend(handles=handles_d, fontsize=9, frameon=True,
            loc='lower center', bbox_to_anchor=(0.5, -0.12), ncols=len(LOC_ORDER))

fig.savefig(os.path.join(FIG_DIR, 'SupplementaryFigure11.png'),
            dpi=600, bbox_inches='tight', facecolor='white')
plt.close()
print("Saved SupplementaryFigure11.png")
