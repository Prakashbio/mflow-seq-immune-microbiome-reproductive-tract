"""Independent numerical/payload checks for the separate v22 figure candidates.

Run after redesign_v22_figures.R. Does not change source measurements or models.
"""
from pathlib import Path
import hashlib, json, sys
import numpy as np
import pandas as pd
from PIL import Image
from pypdf import PdfReader
from lxml import etree

O=Path(sys.argv[1]).resolve()
D=O/'source_data'; A=O/'audit'
checks=[]
def check(name, condition, detail=''):
    checks.append(dict(check=name, passed=bool(condition), detail=detail))
    if not condition: raise AssertionError(name+': '+str(detail))
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
sources=json.loads((A/'Immutable_source_inputs.json').read_text())
check('All staged CSV/RDS inputs retain their original SHA-256',
      all(sha(O/x['file'])==x['sha256'] for x in sources),str(len(sources))+' files')
check('No additional or omitted staged source inputs',
      set(x['file'] for x in sources)==set(str(p.relative_to(O)) for p in D.rglob('*') if p.is_file()))

r=pd.read_csv(D/'FACS_reconstruction/data/Rebuilt_primary_42_woman_site_records.csv')
points=pd.read_csv(A/'Main_Figure_3_points.csv')
s=pd.read_csv(A/'Main_Figure_3_intervals.csv')
c=pd.read_csv(A/'Main_Figure_3_paired_contrasts.csv')
check('FACS primary woman/site records unique',len(r)==42 and not r.duplicated(['Participant','Location']).any())
eps=['IgA_only','IgM_only','Double_positive','IgG_positive']
for ep in eps:
    check(ep+' percentages reproduce recorded event-count denominator',
          np.allclose(r[ep+'_percent'],100*r[ep+'_events']/r.DAPI_negative_events,rtol=0,atol=1e-12))
    for site,g in r.groupby('Location'):
        q=s[(s.Location==site)&(s.endpoint==ep)].iloc[0]
        check(site+' '+ep+' summary',int(q.n_women)==len(g) and
              np.isclose(q.mean_percent,g[ep+'_percent'].mean(),rtol=0,atol=1e-12) and
              np.isclose(q.sample_SD_percent,g[ep+'_percent'].std(ddof=1),rtol=0,atol=1e-12))
    z=points[points.endpoint==ep].merge(r,on=['Participant','Location'],validate='one_to_one')
    check(ep+' displayed points',len(z)==42 and np.allclose(z.percent,z[ep+'_percent'],rtol=0,atol=1e-12))
for site,g in r.groupby('Location'):
    q=c[c.Location.str.upper()==site].iloc[0]
    check(site+' paired contrast',np.isclose(q.mean_M_minus_A_pp,
          (g.IgM_only_percent-g.IgA_only_percent).mean(),rtol=0,atol=1e-12))
original_intervals=pd.read_csv(D/'FACS_reconstruction/data/Rebuilt_Table_2_and_pointwise_CIs.csv')
check('FACS saved CIs unchanged',np.allclose(s[['ci95_low_percent','ci95_high_percent']],
       original_intervals[['ci95_low_percent','ci95_high_percent']],rtol=0,atol=1e-12))
original_c=pd.read_csv(D/'FACS_reconstruction/data/FACS_paired_IgM_minus_IgA_v12.csv')
check('FACS paired CIs unchanged',np.allclose(c[['ci95_low_pp','ci95_high_pp']],
       original_c[original_c['subset']=='All'][['ci95_low_pp','ci95_high_pp']],rtol=0,atol=1e-12))
check('FACS 0–50% display does not hide observations',points.percent.between(0,50).all())

h=pd.read_csv(A/'Main_Figure_6_abundance_panel.csv')
scores=pd.read_csv(D/'taxon_statistics/Site_specific_matched_IgA_IgM_representation_scores.csv')
keys=h.Taxon.unique()
z=scores[(scores.taxon_arm=='baseline_39')&(scores['rank']=='provisional_species')&
         (scores.scope=='All')&scores.Taxon.isin(keys)].copy()
check('Taxon heatmap complete 8 × 3 × 2 cells',len(h)==48 and len(z)==24 and not h.duplicated(['Taxon','Location','Fraction']).any())
for frac in ['IgA','IgM']:
    t=h[h.Fraction==frac].copy();t['Location']=t.Location.str.upper()
    j=t.merge(z,on=['Taxon','Location'],validate='one_to_one')
    check(frac+' heatmap saved means',len(j)==24 and np.allclose(j.pct,100*j['mean_'+frac+'_relative_abundance'],rtol=0,atol=1e-11))
check('Main 6 ratio matches displayed corresponding means',np.allclose(
      z.Score_ratio_of_matched_woman_means,z.mean_IgA_relative_abundance/
      (z.mean_IgA_relative_abundance+z.mean_IgM_relative_abundance),rtol=0,atol=1e-12,equal_nan=True))

m=pd.read_csv(D/'storage_statistics/Storage_library_metrics.csv');m=m[m.taxon_arm=='baseline_39']
v=pd.read_csv(A/'Optional_S18_plotted_values.csv')
check('Optional alpha display sample structure',len(m)==149 and m.Participant.nunique()==15 and m.SpecimenKey.nunique()==40 and len(v)==596)
metricmap={'Shannon entropy':'Shannon','Simpson diversity (1-D)':'Simpson_1_minus_D',
           'Observed richness':'Richness','Pielou evenness':'Pielou_evenness'}
for label,col in metricmap.items():
    j=v[v.metric==label].merge(m,on='SourceCode',validate='one_to_one',suffixes=('_plot','_source'))
    check(label+' values',len(j)==149 and np.allclose(j.value,j[col],rtol=0,atol=1e-12,equal_nan=True))
check('Undefined evenness retained, not replaced by zero',v.value.isna().sum()==6 and
      v[v.value.isna()].metric.eq('Pielou evenness').all() and
      m[m.Pielou_evenness.isna()].Richness.eq(1).all())
v[v.value.isna()].to_csv(A/'Undefined_evenness_records.csv',index=False)
g=scores[(scores.taxon_arm=='baseline_39')&(scores['rank']=='genus_aggregation')&(scores.scope=='All')]
undefined=g[g.Score_ratio_of_matched_woman_means.isna()]
check('Undefined genus scores are zero/zero means',len(undefined)==10 and
      undefined.mean_IgA_relative_abundance.eq(0).all() and undefined.mean_IgM_relative_abundance.eq(0).all())
undefined.to_csv(A/'Undefined_genus_score_records.csv',index=False)
own=pd.read_csv(D/'statistics/baseline_39_Matched_own_PRE.csv')
check('Own-PRE panel denominators',len(own)==104 and own.SpecimenKey.nunique()==37 and own.Participant.nunique()==15)
theme=pd.read_csv(A/'Theme_only_coordinate_preservation.csv')
check('R layer-coordinate comparisons including composite components',len(theme)==20 and theme['pass'].all())

pdfs=list((O/'figures').glob('*.pdf'))
check('20 regenerated figures have PDF, SVG and 600-dpi PNG',len(pdfs)==20 and
      len(list((O/'figures').glob('*.svg')))==20 and
      all(p.with_suffix('.svg').exists() and p.with_suffix('.png').exists() for p in pdfs))
for p in pdfs:
    reader=PdfReader(p)
    check(p.stem+' single-page PDF',len(reader.pages)==1)
    etree.parse(str(p.with_suffix('.svg')))
    with Image.open(p.with_suffix('.png')) as im:
        check(p.stem+' raster resolution',im.width>=4200 and min(im.info.get('dpi',(0,0)))>=599)

# When run inside the full project, verify the sealed v21 payloads and archives.
R=O.parent
if (R/'Simple_submission_v21_2026-11-10/SHA256SUMS.json').exists():
    base=R/'Simple_submission_v21_2026-11-10'
    sealed=json.loads((base/'SHA256SUMS.json').read_text())
    check('All sealed v21 payloads remain unchanged',all(sha(base/k)==val for k,val in sealed.items()),str(len(sealed))+' payloads')
    for name in ['Main_Figure_1.png','Supplementary_Figure_S1.emf','Supplementary_Figure_S2.emf','Supplementary_Figure_S3.emf']:
        check(name+' unchanged asset',sha(O/'figures'/name)==sha(base/'03_figures'/name))
    receipt=json.loads((R/'V21_complete_handoff_receipt_2026-11-10.json').read_text())
    check('Complete v21/v20.1 archive unchanged',sha(R/receipt['archive'])==receipt['sha256'])
    for item in receipt['source_archives']:
        check(item['source_archive']+' unchanged',sha(R/item['source_archive'])==item['source_sha256'])
pd.DataFrame(checks).to_csv(A/'Independent_numerical_and_payload_checks.csv',index=False)
(A/'Audit_summary.json').write_text(json.dumps(dict(status='PASS',checks=len(checks),
    source_files=len(sources),regenerated_figures=20,unchanged_figure_assets=4,
    new_statistical_tests=0,new_fitted_models=0,undefined_evenness=6,undefined_genus_scores=10,
    limitations=['Figure design candidates only; not integrated into Word/PDF submission masters.',
    'Historical gate percentages are not a reconstruction from verified FACSDiva workspace coordinates.',
    'No new experimental DAPI concordance, antibody-specificity or post-sort-purity evidence.']),indent=2)+'\n')
print('PASS',len(checks),'independent checks')
