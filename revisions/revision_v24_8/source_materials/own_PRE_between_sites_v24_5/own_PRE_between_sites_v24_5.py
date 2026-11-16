"""Exploratory paired between-site own-PRE contrasts. No source files are modified.
Run: python own_PRE_between_sites_v24_5.py --inputs INPUT_FOLDER --output OUTPUT_FOLDER
INPUT_FOLDER: baseline_39_source_counts.csv, sensitivity_41_source_counts.csv,
source_metadata_original_codes.csv, prior_within_site_contrasts.csv.
"""
from pathlib import Path
import argparse,hashlib,json,platform,itertools
import numpy as np
import pandas as pd

def bh(p):
 p=np.asarray(p,float);i=np.argsort(p);q=np.empty(len(p));q[i]=np.minimum.accumulate((p[i]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1].clip(0,1);return q

def run(inp,out):
 out.mkdir(parents=True,exist_ok=True);m=pd.read_csv(inp/'source_metadata_original_codes.csv',dtype={'Participant':str});prior=pd.read_csv(inp/'prior_within_site_contrasts.csv');res=[];manifest=[];agg=[];distances=[];draws={};specs={};seed=20261112
 for arm in ['baseline_39','sensitivity_41']:
  x=pd.read_csv(inp/(arm+'_source_counts.csv'),index_col=0);assert list(x.columns)==list(m.SourceCode);assert (x.values>=0).all() and (x.sum()>0).all();rel=x/x.sum(axis=0)
  rows=[]
  for z in m.itertuples():
   if z.Fraction=='PRE':continue
   pre=m[(m.SpecimenKey==z.SpecimenKey)&(m.Fraction=='PRE')];assert len(pre)<=1
   if not len(pre):continue
   pc=pre.iloc[0].SourceCode;a=rel[z.SourceCode].to_numpy();b=rel[pc].to_numpy();d=float(np.abs(a-b).sum()/(a+b).sum());assert 0<=d<=1
   rows.append(dict(taxon_arm=arm,SourceCode=z.SourceCode,PRE=pc,Participant=z.Participant,SpecimenKey=z.SpecimenKey,Location=z.Location,Fraction=z.Fraction,Phase=z.Phase,Distance=d))
  data=pd.DataFrame(rows);assert len(data)==104 and data.Participant.nunique()==15 and data.SpecimenKey.nunique()==37;distances+=rows
  for scope in ['All','Proliferative']:
   z=data if scope=='All' else data[data.Phase=='Proliferative']
   u=z.groupby(['Participant','Location','Fraction'],sort=True).agg(Distance=('Distance','mean'),preparations=('SpecimenKey','nunique'),SourceCodes=('SourceCode',lambda s:';'.join(sorted(s))),PREs=('PRE',lambda s:';'.join(sorted(s)))).reset_index()
   for rr in u.to_dict('records'):agg.append(dict(taxon_arm=arm,scope=scope,**rr))
   rr_scope=[]
   for frac in ['IgA','IgM','IgG']:
    for site1,site2 in [('ENDOMETRIUM','CERVIX'),('ENDOMETRIUM','VAGINA'),('CERVIX','VAGINA')]:
     a=u[(u.Fraction==frac)&(u.Location==site1)];b=u[(u.Fraction==frac)&(u.Location==site2)];j=a.merge(b,on=['Participant','Fraction'],suffixes=('_site1','_site2')).sort_values('Participant',key=lambda s:s.astype(int));v=(j.Distance_site1-j.Distance_site2).to_numpy();n=len(v);assert n>=1
     key=f'{arm}__{scope}__{frac}__{site1}_minus_{site2}';signs=np.array(list(itertools.product([-1,1],repeat=n)),dtype=np.int8);null=signs@v/n;observed=float(v.mean());tol=1e-12;p=float(np.mean(np.abs(null)>=abs(observed)-tol)) if n>=2 else float('nan');s=seed+len(res)+len(rr_scope);rng=np.random.default_rng(s);ix=rng.integers(0,n,size=(5000,n));boot=v[ix].mean(axis=1);lo,hi=np.quantile(boot,[.025,.975],method='linear') if n>=2 else [float('nan'),float('nan')];draws[key+'__indices']=ix;draws[key+'__means']=boot;draws[key+'__signs']=signs;draws[key+'__null_means']=null
     for t in j.to_dict('records'):manifest.append(dict(comparison=key,taxon_arm=arm,scope=scope,site1=site1,site2=site2,**t,difference=t['Distance_site1']-t['Distance_site2']))
     rr_scope.append(dict(comparison=key,taxon_arm=arm,scope=scope,Fraction=frac,site1=site1,site2=site2,women=n,women_ids=';'.join(j.Participant),preparations_site1=int(j.preparations_site1.sum()),preparations_site2=int(j.preparations_site2.sum()),mean_difference=observed,CI95_low=float(lo),CI95_high=float(hi),p_exact_signflip=p,sign_assignments=len(signs),bootstrap_draws=5000,seed=s,inference_status='Exploratory' if n>=2 else 'Descriptive only: one paired woman; no interval or test'))
   newp=np.array([r['p_exact_signflip'] for r in rr_scope]);q=bh(np.nan_to_num(newp,nan=1.));old=prior[(prior.taxon_arm==arm)&(prior.scope==scope)].p_exact_woman_signflip.to_numpy();assert len(old)==9
   combined=bh(list(old)+list(np.nan_to_num(newp,nan=1.)))[9:]
   for r,q9,q18 in zip(rr_scope,q,combined):r.update(q_BH_9_between_site=float(q9) if r['women']>=2 else float('nan'),q_BH_18_all_own_PRE=float(q18) if r['women']>=2 else float('nan'));res.append(r)
 pd.DataFrame(res).to_csv(out/'Between_site_own_PRE_contrasts.csv',index=False);pd.DataFrame(manifest).to_csv(out/'Paired_woman_site_manifest.csv',index=False);pd.DataFrame(agg).to_csv(out/'Woman_site_fraction_averages.csv',index=False);pd.DataFrame(distances).to_csv(out/'Recomputed_own_PRE_distances.csv',index=False);np.savez_compressed(out/'Saved_bootstrap_and_exact_signs.npz',**draws)
 p=pd.DataFrame(res).query("taxon_arm=='baseline_39' and scope=='All'")
 info=dict(analysis='Exploratory between-site difference in distance of the same fraction to its own local PRE',primary='39-assignment full cohort; nine paired site contrasts',unit='Woman; average repeated preparations within woman/site/fraction before matching',estimand='Equal-woman mean of site1 own-PRE distance minus site2 own-PRE distance among women with both sites for that fraction',test='Two-sided complete sign enumeration including observed assignment; assumes independently sampled women and sign-symmetric/exchangeable paired differences under null. Not randomized-site causal inference.',uncertainty='Pointwise percentile woman-bootstrap intervals, 5000 draws per contrast; NumPy linear quantile. Intervals are not simultaneous and do not replace BH.',multiplicity='BH nine new between-site contrasts within each taxon-arm/cohort scope; also BH18 over nine historical within-site fraction and nine between-site tests. Unestimable single-woman comparisons enter adjustment as p=1, but no inferential p/q/CI is reported for them. Sensitivities are not independent replication.',limitations='Missing pairs not imputed; eligible women differ by fraction/site contrast. Anatomy is observational; storage, processing and missingness may confound site-associated differences. Source taxonomy/metadata limitations persist. No inference of antibody-driven effects, viability or storage equivalence.',source_hashes={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in sorted(inp.glob('*.csv'))},software=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__),contrasts=len(res),primary_n_range=[int(p.women.min()),int(p.women.max())],primary_p_range=[float(p.p_exact_signflip.min()),float(p.p_exact_signflip.max())],primary_q_range=[float(p.q_BH_9_between_site.min()),float(p.q_BH_9_between_site.max())])
 (out/'Analysis_specification.json').write_text(json.dumps(info,indent=2)+'\n');print(p[['Fraction','site1','site2','women','mean_difference','CI95_low','CI95_high','p_exact_signflip','q_BH_9_between_site','q_BH_18_all_own_PRE']].to_string(index=False))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();run(a.inputs,a.output)
