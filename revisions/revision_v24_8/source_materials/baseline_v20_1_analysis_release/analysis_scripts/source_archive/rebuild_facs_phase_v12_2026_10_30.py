"""Exploratory phase sensitivities from PDF-count reconstruction and clinical metadata.

No Kasun statistical output is imported. Phase labels are participant metadata
attributed to Merli; independent visit/preparation correspondence is not proven.
"""
from pathlib import Path
from itertools import combinations,product
import argparse,csv,json,hashlib
import numpy as np
parser=argparse.ArgumentParser();parser.add_argument('--package-root',type=Path);a=parser.parse_args()
sd=Path(__file__).resolve().parent
P=a.package_root or (sd.parent if sd.name=='scripts' else sd.parent/'facs_integrated_working_v12_2026-10-30/FACS_reconstruction')
D=P/'data';I=P/'source_inputs';rng=np.random.default_rng(20261030);B=10000
def rows(p):return list(csv.DictReader(p.open()))
def save(name,rr):
 with (D/name).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rr[0]));w.writeheader();w.writerows(rr)
source=rows(D/'Rebuilt_primary_42_woman_site_records.csv');phase=rows(I/'Participant_phase_from_supplied_report.csv')
lookup={r['Participant']:r for r in phase};assert len(lookup)==16
for r in source:
 r.update(Phase=lookup[r['Participant']]['Phase'],reported_phase=lookup[r['Participant']]['reported_phase']);r['M_minus_A_pp']=float(r['IgM_only_percent'])-float(r['IgA_only_percent'])
assert len(source)==42 and len({(r['Participant'],r['Location']) for r in source})==42
source.sort(key=lambda r:(r['Location'],int(r['Participant'])))
save('FACS_phase_linked_42_records_v12.csv',source)
def bh(rr,key,out):
 ps=np.array([r[key] for r in rr]);order=np.argsort(ps);q=np.empty(len(ps));q[order]=np.minimum(1,np.minimum.accumulate((ps[order]*len(ps)/np.arange(1,len(ps)+1))[::-1])[::-1])
 for r,v in zip(rr,q):r[out]=float(v)
paired=[];des=[]
for subset in ['All','Proliferative']:
 family=[]
 for site in ['VAGINA','CERVIX','ENDOMETRIUM']:
  z=[r for r in source if r['Location']==site and (subset=='All' or r['Phase']==subset)]
  x=np.array([r['M_minus_A_pp'] for r in z]);n=len(x);assert 2<=n<=16
  boot=x[rng.integers(0,n,size=(B,n))].mean(axis=1)
  signs=2*((np.arange(2**n,dtype=np.uint64)[:,None]>>np.arange(n,dtype=np.uint64))&1).astype(float)-1
  sims=signs@x/n
  family.append(dict(subset=subset,Location=site,n_women=n,mean_M_minus_A_pp=float(x.mean()),ci95_low_pp=float(np.quantile(boot,.025)),ci95_high_pp=float(np.quantile(boot,.975)),n_positive=int((x>0).sum()),n_negative=int((x<0).sum()),p_exact_signflip=float(np.mean(abs(sims)>=abs(x.mean())-1e-12)),sign_configurations=len(sims),test_statistic='Mean within-woman IgM-gate minus IgA-gate percentage-point difference',null_assumption='Sign exchangeability/symmetry under the null',CI='Pointwise woman-bootstrap percentile; 10000 draws, seed20261030'))
 bh(family,'p_exact_signflip','q_BH_three_sites_within_subset');paired+=family
save('FACS_paired_IgM_minus_IgA_v12.csv',paired)
comparisons=[]
for site in ['VAGINA','CERVIX','ENDOMETRIUM']:
 z=[r for r in source if r['Location']==site]
 x=np.array([r['M_minus_A_pp'] for r in z]);g=np.array([r['Phase']=='Proliferative' for r in z]);a=x[g];b=x[~g]
 boots=a[rng.integers(0,len(a),(B,len(a)))].mean(axis=1)-b[rng.integers(0,len(b),(B,len(b)))].mean(axis=1)
 options=[]
 for storage in sorted({r['Storage'] for r in z}):
  ids=[i for i,r in enumerate(z) if r['Storage']==storage];n=sum(g[ids]);options.append(list(combinations(ids,int(n))))
 sims=[]
 for pick in product(*options):
  ids=[i for block in pick for i in block];mask=np.zeros(len(z),dtype=bool);mask[ids]=True;sims.append(x[mask].mean()-x[~mask].mean())
 effect=a.mean()-b.mean()
 comparisons.append(dict(Location=site,n_proliferative=len(a),n_secretory=len(b),proliferative_mean_M_minus_A_pp=float(a.mean()),secretory_mean_M_minus_A_pp=float(b.mean()),phase_difference_pp=float(effect),ci95_low_pp=float(np.quantile(boots,.025)),ci95_high_pp=float(np.quantile(boots,.975)),p_exact_storage_stratified=float(np.mean(np.abs(sims)>=abs(effect)-1e-12)),label_configurations=len(sims),test_statistic='Proliferative-minus-secretory difference in mean within-woman M-A',permutation='All phase-label allocations preserving phase counts within storage strata',null_assumption='Phase-label exchangeability conditional on recorded storage',CI='Pointwise bootstrap within phase, not storage-adjusted; 10000 draws, seed20261030'))
bh(comparisons,'p_exact_storage_stratified','q_BH_three_sites');save('FACS_exploratory_phase_comparison_v12.csv',comparisons)
for site in ['VAGINA','CERVIX','ENDOMETRIUM']:
 for ph in ['Proliferative','Secretory']:
  z=[r for r in source if r['Location']==site and r['Phase']==ph]
  for ep in ['IgA_only','IgM_only','Double_positive','IgG_positive']:
   x=np.array([float(r[ep+'_percent']) for r in z]);des.append(dict(Location=site,Phase=ph,endpoint=ep,n_women=len(z),mean_percent=float(x.mean()),sample_SD_percent=float(x.std(ddof=1))))
save('FACS_gate_percentages_by_reported_phase_v12.csv',des)
result=dict(source_woman_site_records=42,clinical_women=16,clinical_phase_women={p:len({r['Participant'] for r in source if r['Phase']==p}) for p in ['Proliferative','Secretory']},paired_tests=6,between_phase_tests=3,bootstrap_draws=B,seed=20261030,clinical_phase_assignment_method_independently_verified=False,preparation_to_visit_correspondence_independently_verified=False,raw_FCS_gates_reconstructed=False,metadata_SHA256=hashlib.sha256((I/'Participant_phase_from_supplied_report.csv').read_bytes()).hexdigest(),interpretation='Exploratory, post-hoc and conditional on recorded phase/storage and original report gates. No hormonal-independence, equivalence or causality claim.')
(D/'FACS_phase_analysis_checks_v12.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
