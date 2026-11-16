"""Descriptive PRE correlations, exact participant accumulation and count QC.

Uses Prakash's archived original-code counts, never Kasun's revised outputs.
No inferential correlation/network, adequacy or FACS validation is performed.
"""
from pathlib import Path
from itertools import combinations
from math import comb
import argparse,hashlib,json,shutil,sys,platform
import numpy as np
import pandas as pd

R=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path);p.add_argument('--output',type=Path,default=R/'supplement_audit_working_v6_2026-10-29');a=p.parse_args();O=a.output
I=O/'supplement_audit_inputs';D=O/'supplement_audit_statistics';I.mkdir(parents=True,exist_ok=True);D.mkdir(exist_ok=True)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
if a.inputs:
 for q in a.inputs.iterdir():
  if q.is_file() and q.resolve()!=(I/q.name).resolve():shutil.copy2(q,I/q.name)
else:
 B=R/'source_filter_integrated_working_v5_2026-10-28/source_filter_inputs'
 for n in ['Counts_39_original_codes.csv','Counts_41_original_codes.csv','Counts_268_original_codes.csv','Metadata_original_codes.csv','Workbook_link_manifest.csv']:shutil.copy2(B/n,I/n)
 S=R/'02_evidence/Public_repository_access_2026-10-26/GitHub_snapshot/scripts'
 for n in ['sfig6_cooccurrence_matrix.py','sfig8_taxon_accumulation.py','sfig11_patient_paired_readdepth.py','shared_utils.py']:shutil.copy2(S/n,I/('Original_'+n))
 (I/'Input_origins.json').write_text(json.dumps(dict(counts_source=str(B),scripts_source=str(S),originals_unchanged=True,analysis_owner='Prakash original workflow; not Kasun revised analysis',mapping_limit='Conditional on named-source/workbook correspondence. Historical derived-heading and upstream identities still require author adjudication.'),indent=2)+'\n')
before={q.name:sha(q) for q in I.iterdir() if q.is_file()}
m=pd.read_csv(I/'Metadata_original_codes.csv',dtype=str).fillna('');assert len(m)==149 and m.SourceCode.is_unique and m.Participant.nunique()==15 and m.SpecimenKey.nunique()==40
sites=['VAGINA','CERVIX','ENDOMETRIUM'];arms=['39','41'];librows=[];profiles=[];persons=[];corrows=[];accum=[];incidence=[];distances=[];denoms=[]
for arm in arms:
 x=pd.read_csv(I/f'Counts_{arm}_original_codes.csv',index_col=0);assert x.columns.tolist()==m.SourceCode.tolist() and len(x)==int(arm) and (x.sum()>0).all()
 ra=x.div(x.sum(),axis=1)
 for row in m.itertuples():
  v=ra[row.SourceCode];pos=v[v>0]
  librows.append(dict(arm=arm,SourceCode=row.SourceCode,Participant=row.Participant,Location=row.Location,Fraction=row.Fraction,SpecimenKey=row.SpecimenKey,retained_assignment_counts=int(x[row.SourceCode].sum()),Shannon=float(-(pos*np.log(pos)).sum()),observed_assignments=len(pos)))
 pre=m[m.Fraction=='PRE'];wp={}
 for (site,person),g in pre.groupby(['Location','Participant']):
  # Equal preparation weight after normalizing each library. Never pool counts.
  v=ra[g.SourceCode].mean(axis=1);wp[site,person]=v
  h=[librows[-len(m)+list(m.SourceCode).index(c)]['Shannon'] for c in g.SourceCode]
  persons.append(dict(arm=arm,Location=site,Participant=person,n_PRE_libraries=len(g),SourceCodes=';'.join(g.SourceCode),mean_preparation_Shannon=float(np.mean(h))))
  for t,val in v.items():profiles.append(dict(arm=arm,Location=site,Participant=person,Taxon=t,mean_relative_abundance=float(val),n_PRE_libraries=len(g),SourceCodes=';'.join(g.SourceCode)))
 for site in sites:
  people=sorted(p for s,p in wp if s==site);mat=pd.DataFrame([wp[site,p] for p in people],index=people)
  n=len(mat);assert n==pre[pre.Location==site].Participant.nunique()
  denoms.append(dict(arm=arm,Location=site,n_PRE_women=n,n_PRE_libraries=len(pre[pre.Location==site]),observed_assignments=int((mat>0).any().sum())))
  for t in mat:
   freq=int((mat[t]>0).sum());incidence.append(dict(arm=arm,Location=site,Taxon=t,n_women_present=freq,n_women=n))
  for t1 in mat:
   for t2 in mat:
    valid=mat[t1].nunique()>1 and mat[t2].nunique()>1
    corrows.append(dict(arm=arm,Location=site,Taxon1=t1,Taxon2=t2,n_women=n,rho=float(np.corrcoef(mat[t1].rank(method='average'),mat[t2].rank(method='average'))[0,1]) if valid else np.nan,defined=valid))
  # Exact expectation over all equally likely k-woman subsets; not an estimate
  # of unseen taxa or read-level rarefaction. No Monte Carlo intervals.
  freq=(mat>0).sum(axis=0)
  for k in range(1,n+1):
   mean=sum(1-(comb(n-int(f),k)/comb(n,k) if n-int(f)>=k else 0) for f in freq if f>0)
   accum.append(dict(arm=arm,Location=site,women_accumulated=k,n_available_women=n,expected_observed_assignments=mean,total_observed_assignments=int((freq>0).sum())))
 for s1,s2 in combinations(sites,2):
  common=sorted(set(p for s,p in wp if s==s1)&set(p for s,p in wp if s==s2))
  for p1 in common:
   for p2 in common:
    u,v=wp[s1,p1],wp[s2,p2]
    distances.append(dict(arm=arm,Site1=s1,Site2=s2,Participant1=p1,Participant2=p2,pair_type='Within woman' if p1==p2 else 'Between women',Bray_Curtis=float(np.abs(u-v).sum()/(u.sum()+v.sum())),n_matched_women=len(common)))
allx=pd.read_csv(I/'Counts_268_original_codes.csv',index_col=0)
libs=pd.DataFrame(librows);libs['exported_268_assignment_counts']=libs.SourceCode.map(allx.sum())
links=pd.read_csv(I/'Workbook_link_manifest.csv',dtype={'SampleID':str});libs=libs.merge(links[['SampleID','workbook_reads_numeric','workbook_reads_cell']].rename(columns={'SampleID':'SourceCode'}),on='SourceCode',validate='many_to_one')
libs['denominator_note']='Retained assignment counts in the stated filtered table; not total FASTQ reads, FACS events or a documented pass/fail QC threshold. Workbook count stage unverified.'
pr=pd.DataFrame(profiles);base=pr[pr.arm=='39'];ranking=base.groupby(['Location','Taxon']).mean_relative_abundance.mean().groupby('Taxon').mean().sort_values(ascending=False,kind='stable')
top=pd.DataFrame(dict(display_ID=range(1,13),Taxon=ranking.index[:12],site_balanced_PRE_mean_relative_abundance=ranking.iloc[:12].values));assert len(top)==12
tables={'Library_assignment_count_metrics':libs,'PRE_equal_woman_profiles':pd.DataFrame(profiles),'PRE_participant_Shannon':pd.DataFrame(persons),'PRE_correlation_all_pairs':pd.DataFrame(corrows),'PRE_incidence':pd.DataFrame(incidence),'PRE_exact_accumulation':pd.DataFrame(accum),'PRE_cross_site_distances':pd.DataFrame(distances),'PRE_group_denominators':pd.DataFrame(denoms),'S6_display_taxon_key':top}
for n,df in tables.items():df.to_csv(D/(n+'.csv'),index=False)
complete=set.intersection(*[set(base.loc[base.Location==s,'Participant']) for s in sites])
summary=dict(status='PASS_DESCRIPTIVE_RECONSTRUCTION',input_hashes=before,libraries=149,women=15,preparations=40,PRE_libraries=37,PRE_women_by_site={s:int(pre[pre.Location==s].Participant.nunique()) for s in sites},PRE_complete_three_site_women=len(complete),complete_women=sorted(complete),repeat_rule='Normalize each PRE library first; average profiles and preparation-level Shannon separately within woman/site. Participant 99 endometrium contributes two PRE libraries but one woman.',S6='Top 12 by equal-site/equal-woman PRE mean, 39 membership; every pair in both memberships supplied as CSV. Undefined constant correlations are missing. No p values, significance borders or interaction/network validation.',S8='Exact expectation over k-woman subsets using incidence in PRE profiles; no interval, unseen richness, read rarefaction or sampling-adequacy claim.',S11='PRE trajectories, retained count distribution, descriptive all-ordered within/between PRE distances and count/Shannon scatter. No FACS events, QC pass line, independent-dyad tests or no-depth-bias validation.',hypothesis_tests=False,limits='Conditional named-source mapping; original gate/sequence/control and source-versus-derived-heading issues remain. Species assignments provisional. Composition-dependent correlations are not ecological interactions; selected/filtered taxon accumulation cannot establish total community coverage.')
assert all(sha(I/n)==h for n,h in before.items())
(O/'Supplement_audit_analysis_provenance.json').write_text(json.dumps(summary,indent=2)+'\n')
(D/'Python_environment.json').write_text(json.dumps(dict(python=sys.version,platform=platform.platform(),numpy=np.__version__,pandas=pd.__version__,rank_correlation='Pearson correlation of average ranks; no p values'),indent=2)+'\n')
print(json.dumps({k:v for k,v in summary.items() if k!='input_hashes'},indent=2))
