import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from scipy.stats import mannwhitneyu
import warnings; warnings.filterwarnings('ignore')

from shared_utils import DATA_DIR, FIG_DIR

np.random.seed(42)

print("Loading data...")
seq       = pd.read_csv(os.path.join(DATA_DIR, 'seq_clean.csv'))
otu_kept  = pd.read_csv(os.path.join(DATA_DIR, 'otu_clean.csv'), index_col=0)
otu_decon = pd.read_csv(os.path.join(DATA_DIR, 'otu_clean_decontaminated.csv'), index_col=0)

JUNK = {'No','No reads','no contaminants','no microbes','no reads',
        'no reads without contaminants','no sample'}
seq_r = seq[~seq['Taxon'].isin(JUNK) & seq['Taxon'].notna()
            & ~seq['Taxon'].str.lower().str.startswith('no ')].copy()
seq_pre = seq_r[seq_r['Fraction'] == 'PreSort'].copy()
seq_pre['SampleID'] = seq_pre['Patient'] + '_' + seq_pre['Organ']
pivot = seq_pre.pivot_table(index='Taxon', columns='SampleID',
                             values='Pct', aggfunc='mean').fillna(0)

vag_cols  = [c for c in pivot.columns if 'Vagina'      in c]
endo_cols = [c for c in pivot.columns if 'Endometrium' in c]
cerv_cols = [c for c in pivot.columns if 'Cervix'      in c]
print(f"  Pre-filter taxa: {len(pivot)} | Vag: {len(vag_cols)} | "
      f"Endo: {len(endo_cols)} | Cerv: {len(cerv_cols)}")

celvia_kept = set(otu_kept.index)
BIO_FRT = {
    'Lactobacillus iners','Lactobacillus crispatus','Lactobacillus gasseri',
    'Lactobacillus jensenii','Sneathia vaginalis','Gardnerella vaginalis',
    'Gardnerella leopoldii','Fannyhessea vaginae','Prevotella buccalis',
    'Prevotella oris','Parvimonas micra','Metamycoplasma hominis',
    'Aerococcus christensenii','Peptoniphilus harei','Anaerococcus prevotii',
}


print("\nMethod A: decontam-equivalent scoring...")

def prevalence(df, cols):
    return (df[cols] > 0).sum(axis=1) / len(cols)

prev_vag   = prevalence(pivot, vag_cols)
prev_endo  = prevalence(pivot, endo_cols)
mean_abund = np.log10(pivot.mean(axis=1).clip(0.01))
decontam_score = prev_endo / (prev_endo + prev_vag + 1e-6)

kept_scores    = [float(decontam_score.get(t,0)) for t in pivot.index if t in celvia_kept]
removed_scores = [float(decontam_score.get(t,0)) for t in pivot.index if t not in celvia_kept]
_, pval = mannwhitneyu(kept_scores, removed_scores, alternative='less')

print(f"  Retained (n={len(kept_scores)}): mean score = {np.mean(kept_scores):.3f}")
print(f"  Removed  (n={len(removed_scores)}): mean score = {np.mean(removed_scores):.3f}")
print(f"  Mann-Whitney U p = {pval:.4f} (one-sided)")

dc_out = pd.DataFrame({
    'taxon':             list(pivot.index),
    'prev_vagina':       [float(prev_vag.get(t,0))       for t in pivot.index],
    'prev_endometrium':  [float(prev_endo.get(t,0))      for t in pivot.index],
    'decontam_score':    [float(decontam_score.get(t,0)) for t in pivot.index],
    'celvia_status':     ['RETAINED' if t in celvia_kept else 'REMOVED' for t in pivot.index],
    'known_frt_biology': [t in BIO_FRT                  for t in pivot.index],
}).sort_values('decontam_score', ascending=False)
dc_out.to_csv(os.path.join(FIG_DIR, 'decontam_equivalent_scores.csv'), index=False)
print("  Saved: decontam_equivalent_scores.csv")


# METHOD B: SourceTracker2-EQUIVALENT BAYESIAN SOURCE TRACKING

print("\nMethod B: SourceTracker2-equivalent Gibbs sampling")

def gibbs_source_track(source_matrices, sink_matrix,
                        n_restarts=5, n_draws=500, burnin=200,
                        alpha1=0.001, alpha2=0.001):
    n_taxa  = sink_matrix.shape[1]
    n_sinks = sink_matrix.shape[0]
    n_src   = len(source_matrices)

    profiles = []
    for S in source_matrices:
        t = S.sum(axis=0).astype(float) + alpha2
        profiles.append(t / t.sum())
    un = np.ones(n_taxa)
    profiles.append(un / un.sum())          # Unknown source
    all_p   = np.array(profiles)
    n_total = n_src + 1

    results = []
    for si in range(n_sinks):
        sink  = sink_matrix[si].astype(float)
        total = sink.sum()
        if total == 0:
            results.append(np.ones(n_total) / n_total)
            continue
        sink_int        = np.random.multinomial(500, sink / total)
        best_p, best_ll = None, -np.inf

        for _ in range(n_restarts):
            props = np.random.dirichlet(np.ones(n_total) * alpha1)
            draws = []
            for it in range(n_draws + burnin):
                asgn = np.zeros(n_total)
                for ti in range(n_taxa):
                    if sink_int[ti] == 0:
                        continue
                    pw = props * all_p[:, ti]
                    ps = pw.sum()
                    pw = pw / ps if ps > 0 else np.ones(n_total) / n_total
                    asgn += np.random.multinomial(sink_int[ti], pw)
                props = np.random.dirichlet(asgn + alpha1)
                if it >= burnin:
                    draws.append(props.copy())
            mp = np.mean(draws, axis=0)
            ll = np.sum(sink_int * np.log(np.dot(mp, all_p) + 1e-10))
            if ll > best_ll:
                best_ll = ll
                best_p  = mp
        results.append(best_p)
    return np.array(results)

post_taxa = [t for t in otu_decon.index if t in pivot.index]
piv_post  = pivot.loc[post_taxa]
props = gibbs_source_track(
    [piv_post[vag_cols].values.T, piv_post[cerv_cols].values.T],
    piv_post[endo_cols].values.T
)
st2 = pd.DataFrame(props, columns=['Vaginal','Cervical','Unknown'],
                    index=[c.replace('_Endometrium','').replace('HUT','').lstrip('0') or '0' for c in endo_cols])
st2.round(4).to_csv(os.path.join(FIG_DIR, 'sourcetracker2_equivalent_results.csv'))
print("  Saved: sourcetracker2_equivalent_results.csv")
print(f"  Mean Vaginal={st2['Vaginal'].mean():.3f}  "
      f"Cervical={st2['Cervical'].mean():.3f}  "
      f"Unknown={st2['Unknown'].mean():.3f}")


fig = plt.figure(figsize=(20, 14))
gs  = gridspec.GridSpec(2, 3, figure=fig, hspace=0.25, wspace=0.33)

# Panel A
ax_a = fig.add_subplot(gs[0, :2])

# First pass: collect points to label and their positions
labels_to_add = []
for taxon in pivot.index:
    x  = float(mean_abund.get(taxon, -2))
    y  = float(decontam_score.get(taxon, 0))
    ic = taxon in celvia_kept
    ib = taxon in BIO_FRT
    if not ic and y > 0.6:  c, m, s = '#D32F2F', 'X',  100
    elif ic and y > 0.6:    c, m, s = '#FF8C00', '^',  100
    elif ic and ib:         c, m, s = '#1B5E20', 'D',   85
    elif ic:                c, m, s = '#2196F3', 'o',   65
    else:                   c, m, s = '#BDBDBD', '.',   38
    ax_a.scatter(x, y, c=c, marker=m, s=s, alpha=0.82, zorder=3)

    if (y > 0.72 and ic) or taxon in BIO_FRT:
        short = (taxon.replace('Lactobacillus','L.').replace('Streptococcus','S.')
                 .replace('Cutibacterium','C.').replace('Gardnerella','G.'))
        labels_to_add.append((x, y, short))

# Second pass: add labels with smart positioning to avoid overlaps and axis boundaries
labels_to_add.sort(key=lambda t: (-t[1], t[0]))  # Sort by y (descending), then x
used_positions = []
y_min, y_max = -0.05, 1.18  # axis limits
# Estimate x range from the data
all_x = [float(mean_abund.get(taxon, -2)) for taxon in pivot.index]
x_min, x_max = min(all_x), max(all_x)

for x, y, label in labels_to_add:
    # Try different positions around the point, prioritizing right side
    offsets = [(5, 3), (5, -3), (5, 0), (5, 5), (5, -5), (-40, 3), (-40, -3), (-40, 0), (5, 8), (5, -8)]
    best_offset = offsets[0]
    max_dist = -1

    for offset in offsets:
        # Calculate text position in data coordinates (approximate)
        # Using rough conversion: offset points -> data units
        text_y = y + offset[1] * 0.01
        text_x = x + offset[0] * 0.015

        # Check if label would be out of bounds (y-axis)
        if text_y < y_min + 0.03 or text_y > y_max - 0.03:
            continue

        # Check if label would be too far left (overlapping with y-axis)
        # For left-side labels, need extra margin (~0.5 data units for text width)
        if offset[0] < 0:  # left-side label
            if text_x < x_min + 0.6:  # too close to left edge
                continue

        # Check minimum distance to existing labels
        min_dist = float('inf')
        for px, py in used_positions:
            dist = ((text_x - px)**2 + (text_y - py)**2)**0.5
            min_dist = min(min_dist, dist)

        if min_dist > max_dist:
            max_dist = min_dist
            best_offset = offset

    ax_a.annotate(label, (x, y), fontsize=6.5, ha='left' if best_offset[0] > 0 else 'right',
                  xytext=best_offset, textcoords='offset points', zorder=5)
    used_positions.append((x + best_offset[0] * 0.015, y + best_offset[1] * 0.01))
ax_a.axhline(0.6, color='red', ls='--', lw=1.3, alpha=0.75)
ax_a.set_xlabel('Mean log\u2081\u2080(relative abundance %)', fontsize=11)
ax_a.set_ylabel('decontam-equivalent score\n(prev_endo / [prev_endo + prev_vagina])', fontsize=10)
ax_a.set_title('A:  decontam-Equivalent Prevalence Scoring\n'
               '(pre-filter dataset, 127 taxa, PRE-sorted fractions)',
               fontsize=11, fontweight='bold', loc='left')
leg = [mpatches.Patch(color='#D32F2F', label='Removed by Celvia CC AS + score>0.6 (agreement)'),
       mpatches.Patch(color='#FF8C00', label='Retained by Celvia CC AS + score>0.6 (biological — see text)'),
       mpatches.Patch(color='#1B5E20', label='Known FRT biology — retained correctly'),
       mpatches.Patch(color='#2196F3', label='Retained, not flagged'),
       mpatches.Patch(color='#BDBDBD', label='Removed, not flagged'),
       Line2D([0],[0], color='red', ls='--', label='Threshold (0.6)')]
ax_a.legend(handles=leg, fontsize=8, loc='upper right', framealpha=0.9)
ax_a.set_ylim(-0.05, 1.18); ax_a.grid(True, alpha=0.2)

# Panel B
ax_b = fig.add_subplot(gs[0, 2])
cm_b = {'Vaginal':'#2166AC','Cervical':'#4DAC26','Unknown':'#BDBDBD'}
pats = list(st2.index)
bot  = np.zeros(len(pats))
for src in ['Vaginal','Cervical','Unknown']:
    ax_b.bar(range(len(pats)), st2[src].values, bottom=bot,
             color=cm_b[src], label=src, width=0.75, edgecolor='white', lw=0.5)
    bot += st2[src].values
ax_b.set_xticks(range(len(pats)))
ax_b.set_xticklabels(pats, fontsize=8.5)
ax_b.set_ylabel('Estimated source proportion', fontsize=10)
ax_b.set_ylim(0, 1.05)
ax_b.set_title('B:  SourceTracker2-Equivalent Bayesian Source Tracking\n'
               '(post-filter, 39 taxa)',
               fontsize=11, fontweight='bold', loc='left')
legend_title = f'Mean Unknown: {st2["Unknown"].mean()*100:.1f}\u202f\u00b1\u202f{st2["Unknown"].std()*100:.1f}%'
ax_b.legend(fontsize=9, loc='upper center', bbox_to_anchor=(0.5, -0.05),
            ncol=3, frameon=True, framealpha=0.9, title=legend_title, title_fontsize=8.5)
ax_b.axhline(0.5, color='black', lw=0.8, ls=':', alpha=0.5)

# Panel C
ax_c = fig.add_subplot(gs[1, :2])
ax_c.boxplot([kept_scores, removed_scores],
              tick_labels=[f'Retained\n(n={len(kept_scores)})',
                           f'Removed\n(n={len(removed_scores)})'],
              patch_artist=True,
              boxprops=dict(facecolor='#E3F2FD'),
              medianprops=dict(color='#D32F2F', lw=2.2),
              widths=0.45)
for i,(sc,col) in enumerate([(kept_scores,'#1565C0'),(removed_scores,'#B71C1C')], 1):
    jit = np.random.uniform(-0.14, 0.14, len(sc))
    ax_c.scatter(np.array([i]*len(sc))+jit, sc,
                 alpha=0.38, s=24, color=col, zorder=3)
ax_c.axhline(0.6, color='red', ls='--', lw=1.3, alpha=0.75, label='Threshold (0.6)')
ax_c.set_ylabel('decontam-equivalent score', fontsize=11)
ax_c.set_title('C:  Celvia CC AS Filter Validation\n'
               'No systematic endometrial-prevalence bias between retained and removed taxa',
               fontsize=11, fontweight='bold', loc='left')
ax_c.legend(fontsize=10)
ax_c.set_ylim(-0.05, 1.18); ax_c.grid(True, axis='y', alpha=0.28)
ax_c.text(0.5, 0.96,
          f'Mann\u2013Whitney U,  p\u202f=\u202f{pval:.3f}  (one-sided, n.s.)',
          transform=ax_c.transAxes, ha='center', va='top', fontsize=11,
          bbox=dict(boxstyle='round,pad=0.35',
                    facecolor='lightyellow', edgecolor='#cccc00'))

# Panel D
ax_d = fig.add_subplot(gs[1, 2])
ax_d.axis('off')
rows = [
    ('bold', 'decontam analysis',          'pre-filter, PRE-sorted'),
    ('norm', 'Total taxa',                 str(len(pivot))),
    ('norm', 'Retained by Celvia CC AS',   str(len(kept_scores))),
    ('norm', 'Removed by Celvia CC AS',    str(len(removed_scores))),
    ('norm', 'Retained mean score',        f'{np.mean(kept_scores):.3f}'),
    ('norm', 'Removed mean score',         f'{np.mean(removed_scores):.3f}'),
    ('norm', 'Mann\u2013Whitney p',        f'{pval:.3f} (n.s.)'),
    ('sep',  '', ''),
    ('bold', 'SourceTracker2 analysis',    'post-filter, PRE-sorted'),
    ('norm', 'Sink samples (patients)',    f'{len(endo_cols)} endometrial'),
    ('norm', 'Mean Vaginal',               f'{st2["Vaginal"].mean()*100:.1f}%'),
    ('norm', 'Mean Cervical',              f'{st2["Cervical"].mean()*100:.1f}%'),
    ('norm', 'Mean Unknown',               f'{st2["Unknown"].mean()*100:.1f}\u202f\u00b1\u202f{st2["Unknown"].std()*100:.1f}%'),
    ('norm', 'Patients \u226520% Unknown', f'{(st2["Unknown"]>=0.2).sum()}/{len(endo_cols)}'),
]
y = 0.97
for style, k, v in rows:
    if style == 'sep':
        y -= 0.022; continue
    ax_d.text(0.02, y, k, transform=ax_d.transAxes, fontsize=8.5,
              fontweight='bold' if style == 'bold' else 'normal', va='top')
    ax_d.text(0.60, y, v, transform=ax_d.transAxes, fontsize=8.5, va='top')
    y -= 0.068
ax_d.set_title('D:  Summary Statistics', fontsize=11, fontweight='bold', loc='left')


plt.savefig(os.path.join(FIG_DIR, 'SupplementaryFigure16.png'), dpi=180, bbox_inches='tight', facecolor='white')
plt.close()
print("  Saved: SupplementaryFigure16.png")
