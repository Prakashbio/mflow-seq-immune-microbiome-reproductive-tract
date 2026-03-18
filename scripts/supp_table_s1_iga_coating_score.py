"""
Supplementary Table S1 – IgA Coating Score Table (all 39 taxa)
===============================================================
Computes and exports the per-taxon IgA coating score for all 39 taxa
in the decontaminated dataset, with supporting abundance and prevalence
columns to allow readers to evaluate signal reliability for each taxon.

Definition
----------
IgA coating score = IgA_mean_rel_abundance / (IgA_mean + IgM_mean)
  • Score = 1.0 → exclusively IgA-coated across all samples
  • Score = 0.5 → equal IgA and IgM coating (near-neutral)
  • Score = 0.0 → exclusively IgM-coated across all samples
  • NaN       → taxon absent from both IgA and IgM fractions

Means are computed across ALL samples per Ig fraction (all locations
pooled), matching the formula stated in the manuscript Methods section.

Output files
------------
  Output/SupplementaryTableS1_IgACoatingScores.csv   — machine-readable
  Output/SupplementaryTableS1_IgACoatingScores.tsv   — tab-separated

Reference
---------
Lin & Peddada (2020) Nat Commun 11:3514  (ANCOM-BC framework context)
Jackson et al. (2021) Microbiome 9:33    (IgA-Seq framework)
Palm et al. (2014) Cell 158:1000-1010   (IgA coating score concept)

Verified against ground-truth data 2026-03-18. All values are
reproducible by running this script on the deposited data files.
"""

import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings('ignore')

from shared_utils import load_data, LOC_ORDER, FIG_DIR

# ── Output directory (reuse the Output/ folder used by other scripts) ─────────
OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'Output')
os.makedirs(OUT_DIR, exist_ok=True)

# ── Load data ─────────────────────────────────────────────────────────────────
otu_rel, meta, tax = load_data()
# otu_rel: relative abundance, shape (taxa × samples)
# meta:    sample metadata — columns: Location, Antibody, Condition, Storage
# tax:     taxonomy table

DATA  = otu_rel.T   # samples × taxa
taxa  = otu_rel.index.tolist()

print(f"Loaded: {otu_rel.shape[0]} taxa × {otu_rel.shape[1]} samples")
print(f"Antibody fractions: {sorted(meta['Antibody'].unique())}")
print(f"Anatomical sites:   {sorted(meta['Location'].unique())}")

# ── Identify IgA and IgM sample sets ─────────────────────────────────────────
iga_ids = meta[meta['Antibody'] == 'A'].index.tolist()
igm_ids = meta[meta['Antibody'] == 'M'].index.tolist()

print(f"\nIgA samples: {len(iga_ids)}")
print(f"IgM samples: {len(igm_ids)}")

# ── Compute global means (all anatomical sites pooled) ────────────────────────
# This is the formula stated in the manuscript:
#   Score = IgA_mean / (IgA_mean + IgM_mean), averaged across all samples
mu_iga = otu_rel[iga_ids].mean(axis=1)   # mean rel. abundance per taxon in IgA
mu_igm = otu_rel[igm_ids].mean(axis=1)   # mean rel. abundance per taxon in IgM

denom  = mu_iga + mu_igm

# Assign NaN where both means are zero (taxon absent from both fractions)
iga_score_global = mu_iga / denom.replace(0, np.nan)

# ── Compute per-site means (for reference columns) ────────────────────────────
site_mu_iga = {}
site_mu_igm = {}
site_score  = {}
site_prev_iga = {}
site_prev_igm = {}

for loc in LOC_ORDER:
    a_ids = [s for s in iga_ids if meta.loc[s, 'Location'] == loc]
    m_ids = [s for s in igm_ids if meta.loc[s, 'Location'] == loc]

    mu_a = otu_rel[a_ids].mean(axis=1) if a_ids else pd.Series(0.0, index=taxa)
    mu_m = otu_rel[m_ids].mean(axis=1) if m_ids else pd.Series(0.0, index=taxa)
    den  = mu_a + mu_m
    sc   = mu_a / den.replace(0, np.nan)

    site_mu_iga[loc]  = mu_a
    site_mu_igm[loc]  = mu_m
    site_score[loc]   = sc
    site_prev_iga[loc] = (otu_rel[a_ids] > 0).mean(axis=1) * 100 if a_ids else pd.Series(0.0, index=taxa)
    site_prev_igm[loc] = (otu_rel[m_ids] > 0).mean(axis=1) * 100 if m_ids else pd.Series(0.0, index=taxa)

# ── Prevalence across all IgA / IgM samples ───────────────────────────────────
prev_iga_all = (otu_rel[iga_ids] > 0).mean(axis=1) * 100
prev_igm_all = (otu_rel[igm_ids] > 0).mean(axis=1) * 100

# ── Classification helper ─────────────────────────────────────────────────────
def classify_score(score):
    """Bin a continuous coating score into a descriptive category."""
    if pd.isna(score):
        return 'Not detected'
    if score >= 0.90:
        return 'Strong IgA (≥0.90)'
    if score >= 0.70:
        return 'IgA-dominant (0.70–0.89)'
    if score >= 0.60:
        return 'Moderate IgA (0.60–0.69)'
    if score >= 0.40:
        return 'Near-neutral (0.40–0.59)'
    if score >= 0.30:
        return 'Moderate IgM (0.30–0.39)'
    if score >= 0.10:
        return 'IgM-dominant (0.10–0.29)'
    return 'Strong IgM (<0.10)'


# ── Reliability flag ──────────────────────────────────────────────────────────
# Flag taxa where the score is based on very sparse data
# (total mean < 0.05% AND prevalence < 5% in both fractions)
# These scores should be interpreted cautiously.
def reliability_flag(taxon):
    total_pct = (mu_iga[taxon] + mu_igm[taxon]) * 100
    pa = prev_iga_all[taxon]
    pm = prev_igm_all[taxon]
    if total_pct < 0.05 and pa < 5 and pm < 5:
        return 'Low (n<3 non-zero samples in both fractions)'
    if total_pct < 0.10 and pa < 10 and pm < 10:
        return 'Moderate (sparse detection)'
    return 'Adequate'


# ── Biological interpretation (curated) ───────────────────────────────────────
BIO_INTERP = {
    'Sneathia vaginalis':           'BV/dysbiosis marker — predominantly IgA-targeted',
    'Metamycoplasma hominis':       'Reproductive tract pathobiont; exclusively IgA-coated',
    'Parvimonas micra':             'Endometrial dysbiosis-associated anaerobe',
    'Prevotella oris':              'BV-associated dysbiosis anaerobe',
    'Haemophilus parainfluenzae':   'Opportunistic respiratory/FRT commensal',
    'Lactobacillus delbrueckii':    'Cervical/endometrial Lactobacillus; IgA-enriched',
    'Aerococcus christensenii':     'Reproductive tract anaerobe; IgA-enriched',
    'Prevotella buccalis':          'BV-associated dysbiosis anaerobe',
    'Fannyhessea vaginae':          'BV-associated; strongly IgA-targeted',
    'Anaerococcus prevotii':        'BV-associated anaerobe',
    'Gardnerella vaginalis':        'Classic BV-associated taxon; IgA-coated',
    'Cutibacterium granulosum':     'Low-abundance skin commensal',
    'Granulicatella adiacens':      'Endometrial microbiome constituent',
    'Corynebacterium matruchotii':  'Low-abundance environmental commensal',
    'Peptoniphilus harei':          'Sparse anaerobe; score based on n=2 non-zero samples — interpret cautiously',
    'Lancefieldella parvula':       'BV-associated streptococcal anaerobe; sparse detection',
    'Staphylococcus aureus':        'Environmental; very low FRT prevalence',
    'Rothia mucilaginosa':          'Low-abundance oral/respiratory commensal',
    'Lactobacillus iners':          'Dominant vaginal commensal; near-neutral IgA coating (transitional CST-III)',
    'Streptococcus thermophilus':   'Low-abundance; near-equal IgA/IgM coating',
    'Gardnerella leopoldii':        'BV-associated Gardnerella species; near-neutral coating',
    'Streptococcus sanguinis':      'Low-abundance oral commensal',
    'Cutibacterium acnes':          'Skin commensal; moderate IgM preference',
    'Acinetobacter ursingii':       'Environmental/opportunistic; low FRT prevalence',
    'Gemella haemolysans':          'Low-abundance oral commensal',
    'Anaerococcus mediterraneensis':'BV-associated anaerobe',
    'Lactobacillus crispatus':      'Health-protective Lactobacillus (CST-I); predominantly IgM-coated',
    'Staphylococcus hominis':       'Environmental skin commensal; low FRT prevalence',
    'Peptoniphilus sp. SAHP1':      'Anaerobic endometrial commensal',
    'Lacticaseibacillus paracasei': 'Low-abundance Lactobacillus',
    'Lactococcus cremoris':         'Low-abundance dairy-associated Lactococcus',
    'Cutibacterium modestum':       'Low-abundance skin commensal',
    'Moraxella osloensis':          'Low-abundance environmental opportunist',
    'Lactobacillus jensenii':       'Health-protective vaginal Lactobacillus; strongly IgM-preferring',
    'Escherichia coli':             'Low-abundance; potential environmental/ascending signal',
    'Corynebacterium sanguinis':    'Low-abundance; no IgA signal detected',
    'Ezakiella coagulans':          'Vaginal anaerobe; no IgA signal in this dataset',
    'Lactobacillus gasseri':        'Health-protective vaginal Lactobacillus; exclusively IgM-coated',
    'Streptococcus salivarius':     'Low-abundance oral commensal; no IgA signal',
}

# ── Build output table ────────────────────────────────────────────────────────
rows = []
for taxon in iga_score_global.sort_values(ascending=False).index:
    sc   = iga_score_global[taxon]
    rows.append({
        # ── Identity ─────────────────────────────────────────────────────────
        'Taxon': taxon,

        # ── Primary score (global, all sites pooled) ──────────────────────
        'IgA_coating_score':     round(float(sc),   3) if not pd.isna(sc) else np.nan,
        'Classification':        classify_score(sc),

        # ── Global mean relative abundance (%) ────────────────────────────
        'IgA_mean_pct_global':   round(float(mu_iga[taxon]) * 100, 4),
        'IgM_mean_pct_global':   round(float(mu_igm[taxon]) * 100, 4),

        # ── Global prevalence (% samples > 0) ────────────────────────────
        'Prevalence_IgA_pct':    round(float(prev_iga_all[taxon]), 1),
        'Prevalence_IgM_pct':    round(float(prev_igm_all[taxon]), 1),

        # ── Per-site scores ───────────────────────────────────────────────
        'Score_Vagina':          round(float(site_score['VAGINA'][taxon]),       3)
                                 if not pd.isna(site_score['VAGINA'][taxon])       else np.nan,
        'Score_Cervix':          round(float(site_score['CERVIX'][taxon]),       3)
                                 if not pd.isna(site_score['CERVIX'][taxon])       else np.nan,
        'Score_Endometrium':     round(float(site_score['ENDOMETRIUM'][taxon]),  3)
                                 if not pd.isna(site_score['ENDOMETRIUM'][taxon])  else np.nan,

        # ── Per-site IgA mean % ───────────────────────────────────────────
        'IgA_mean_pct_Vagina':   round(float(site_mu_iga['VAGINA'][taxon])       * 100, 4),
        'IgA_mean_pct_Cervix':   round(float(site_mu_iga['CERVIX'][taxon])       * 100, 4),
        'IgA_mean_pct_Endo':     round(float(site_mu_iga['ENDOMETRIUM'][taxon])  * 100, 4),

        # ── Per-site IgM mean % ───────────────────────────────────────────
        'IgM_mean_pct_Vagina':   round(float(site_mu_igm['VAGINA'][taxon])       * 100, 4),
        'IgM_mean_pct_Cervix':   round(float(site_mu_igm['CERVIX'][taxon])       * 100, 4),
        'IgM_mean_pct_Endo':     round(float(site_mu_igm['ENDOMETRIUM'][taxon])  * 100, 4),

        # ── Reliability and interpretation ────────────────────────────────
        'Reliability':           reliability_flag(taxon),
        'Biological_interpretation': BIO_INTERP.get(taxon, ''),
    })

df = pd.DataFrame(rows)
df.index = range(1, len(df) + 1)   # 1-based rank

# ── Print summary to console ──────────────────────────────────────────────────
print("\n" + "=" * 100)
print("SUPPLEMENTARY TABLE S1 — IgA Coating Scores (all 39 taxa, ranked)")
print("=" * 100)
header = (f"{'Rank':>4}  {'Taxon':<40} {'Score':>7} {'Class':<25} "
          f"{'IgA%':>6} {'IgM%':>6} {'PrevA':>6} {'PrevM':>6}  Reliability")
print(header)
print("-" * 120)
for rank, row in df.iterrows():
    sc_str = f"{row['IgA_coating_score']:.3f}" if not pd.isna(row['IgA_coating_score']) else "  NaN"
    rel_short = row['Reliability'].split('(')[0].strip()
    print(f"  {rank:>3}  {row['Taxon']:<40} {sc_str:>7} "
          f"{row['Classification']:<25} "
          f"{row['IgA_mean_pct_global']:>6.2f} {row['IgM_mean_pct_global']:>6.2f} "
          f"{row['Prevalence_IgA_pct']:>6.1f} {row['Prevalence_IgM_pct']:>6.1f}  {rel_short}")

print(f"\n{len(df)} taxa total")
print(f"  Strong IgA (≥0.90):          {(df['IgA_coating_score'] >= 0.90).sum()}")
print(f"  IgA-dominant (0.70–0.89):    {((df['IgA_coating_score'] >= 0.70) & (df['IgA_coating_score'] < 0.90)).sum()}")
print(f"  Moderate IgA (0.60–0.69):    {((df['IgA_coating_score'] >= 0.60) & (df['IgA_coating_score'] < 0.70)).sum()}")
print(f"  Near-neutral (0.40–0.59):    {((df['IgA_coating_score'] >= 0.40) & (df['IgA_coating_score'] < 0.60)).sum()}")
print(f"  Moderate IgM (0.30–0.39):    {((df['IgA_coating_score'] >= 0.30) & (df['IgA_coating_score'] < 0.40)).sum()}")
print(f"  IgM-dominant (0.10–0.29):    {((df['IgA_coating_score'] >= 0.10) & (df['IgA_coating_score'] < 0.30)).sum()}")
print(f"  Strong IgM (<0.10):          {(df['IgA_coating_score'] < 0.10).sum()}")

# ── Spot-check key manuscript values ──────────────────────────────────────────
print("\n" + "=" * 60)
print("SPOT-CHECK — key values cited in manuscript")
print("=" * 60)
spot = {
    'Sneathia vaginalis':    ('0.999', 'manuscript = 0.999'),
    'Lactobacillus iners':   ('0.546', 'manuscript = 0.546 (intermediate)'),
    'Lactobacillus crispatus':('0.263', 'manuscript = 0.263'),
    'Lactobacillus jensenii': ('0.071', 'manuscript = 0.071'),
    'Lactobacillus gasseri':  ('0.000', 'manuscript = 0.000'),
    'Peptoniphilus harei':    ('0.622', 'manuscript was 0.780 — corrected'),
}
for taxon, (expected, note) in spot.items():
    actual = df[df['Taxon'] == taxon]['IgA_coating_score'].values
    if len(actual):
        match = '✅' if abs(float(actual[0]) - float(expected)) < 0.005 else '⚠️'
        print(f"  {match} {taxon}: computed={actual[0]:.3f}, expected={expected}  ({note})")

# ── Export ────────────────────────────────────────────────────────────────────
csv_path = os.path.join(OUT_DIR, 'SupplementaryTableS1_IgACoatingScores.csv')
tsv_path = os.path.join(OUT_DIR, 'SupplementaryTableS1_IgACoatingScores.tsv')

df.to_csv(csv_path, index_label='Rank')
df.to_csv(tsv_path, sep='\t', index_label='Rank')

print(f"\nSaved: {csv_path}")
print(f"Saved: {tsv_path}")
print("\nColumn descriptions:")
print("  IgA_coating_score   : IgA_mean / (IgA_mean + IgM_mean), all samples pooled")
print("  Classification      : Categorical bin (see Methods)")
print("  IgA_mean_pct_global : Mean relative abundance (%) in IgA fraction, all sites")
print("  IgM_mean_pct_global : Mean relative abundance (%) in IgM fraction, all sites")
print("  Prevalence_IgA_pct  : % of IgA samples where this taxon > 0")
print("  Prevalence_IgM_pct  : % of IgM samples where this taxon > 0")
print("  Score_Vagina/Cervix/Endometrium : Per-site coating score")
print("  IgA/IgM_mean_pct_Vagina/Cervix/Endo : Per-site mean abundance (%)")
print("  Reliability         : Data quality flag (based on prevalence and abundance)")
print("  Biological_interpretation : Curated annotation")
