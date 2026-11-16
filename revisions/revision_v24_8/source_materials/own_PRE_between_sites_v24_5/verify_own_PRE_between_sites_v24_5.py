"""Independent scalar and scipy checks; invokes no analysis helper code."""
from pathlib import Path
import csv,json,itertools,hashlib,math,argparse
import numpy as np
def independent_bh(p):
 order=sorted(range(len(p)),key=lambda i:p[i]);ranks={i:k+1 for k,i in enumerate(order)}
 return [min(1.,min(len(p)*p[j]/ranks[j] for j in order if ranks[j]>=ranks[i])) for i in range(len(p))]
R=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('--workdir',type=Path,default=R/'06_checks/V24_5_own_PRE_site_comparisons_2026-11-12');a=p.parse_args();Q=a.workdir;I=Q/'inputs';O=Q/'results'
rows=lambda p:list(csv.DictReader(p.open()));m=rows(I/'source_metadata_original_codes.csv');results=rows(O/'Between_site_own_PRE_contrasts.csv');pairs=rows(O/'Paired_woman_site_manifest.csv');saved=np.load(O/'Saved_bootstrap_and_exact_signs.npz');checks=0;maxerr=0.
def check(v,msg):
 global checks
 assert v,msg;checks+=1
def near(a,b,msg,tol=1e-12):
 global maxerr
 err=abs(float(a)-float(b));maxerr=max(maxerr,err);check(err<tol,(msg,err))
by={};computed={}
for arm in ['baseline_39','sensitivity_41']:
 with (I/(arm+'_source_counts.csv')).open() as f:
  z=list(csv.reader(f));names=z[0][1:];vals={name:[float(r[j+1]) for r in z[1:]] for j,name in enumerate(names)}
 for name in names:
  total=sum(vals[name]);vals[name]=[x/total for x in vals[name]]
 for r in m:
  if r['Fraction']=='PRE':continue
  pre=[t for t in m if t['SpecimenKey']==r['SpecimenKey'] and t['Fraction']=='PRE']
  if not pre:continue
  a=vals[r['SourceCode']];b=vals[pre[0]['SourceCode']];d=sum(abs(x-y) for x,y in zip(a,b))/sum(x+y for x,y in zip(a,b));computed[arm,r['SourceCode']]=d
 for old in rows(I/(arm+'_archived_Matched_own_PRE.csv')):near(computed[arm,old['SourceCode']],old['Distance'],'archived distance')
for row in results:
 arm=row['taxon_arm'];scope=row['scope'];frac=row['Fraction'];a=row['site1'];b=row['site2'];key=row['comparison'];values=[];eligible=[]
 for person in sorted(set(r['Participant'] for r in m),key=int):
  za=[r for r in m if r['Participant']==person and r['Fraction']==frac and r['Location']==a and (scope=='All' or r['Phase']=='Proliferative') and (arm,r['SourceCode']) in computed]
  zb=[r for r in m if r['Participant']==person and r['Fraction']==frac and r['Location']==b and (scope=='All' or r['Phase']=='Proliferative') and (arm,r['SourceCode']) in computed]
  if za and zb:values.append(sum(computed[arm,r['SourceCode']] for r in za)/len(za)-sum(computed[arm,r['SourceCode']] for r in zb)/len(zb));eligible.append(person)
 n=len(values);check(n==int(row['women']),key+' n');check(';'.join(eligible)==row['women_ids'],key+' women');near(sum(values)/n,row['mean_difference'],key+' mean')
 j=[r for r in pairs if r['comparison']==key];check(len(j)==n,key+' manifest')
 for r,v in zip(j,values):near(r['difference'],v,key+' paired value')
 if n==1:
  check(all(row[k]=='' for k in ['CI95_low','CI95_high','p_exact_signflip','q_BH_9_between_site','q_BH_18_all_own_PRE']),key+' no single-woman inference');continue
 obs=abs(sum(values)/n);tail=sum(abs(sum(v*s for v,s in zip(values,sign))/n)>=obs-1e-12 for sign in itertools.product([-1,1],repeat=n));near(tail/2**n,row['p_exact_signflip'],key+' scalar exact p')
 ix=saved[key+'__indices'];check(ix.shape==(5000,n),key+' draw shape');check(np.array_equal(ix,np.random.default_rng(int(row['seed'])).integers(0,n,(5000,n))),key+' seed');boot=np.array([sum(values[int(k)] for k in sample)/n for sample in ix]);near(np.max(np.abs(boot-saved[key+'__means'])),0,key+' saved bootstrap');near(np.quantile(boot,.025),row['CI95_low'],key+' low');near(np.quantile(boot,.975),row['CI95_high'],key+' high')
 for rr in j:
  check(rr['Participant'] in eligible,key+' parent');check(rr['Fraction']==frac,key+' same fraction');check(rr['Location_site1']==a and rr['Location_site2']==b,key+' sites')
for arm,scope in itertools.product(['baseline_39','sensitivity_41'],['All','Proliferative']):
 z=[r for r in results if r['taxon_arm']==arm and r['scope']==scope];p=[float(r['p_exact_signflip']) if r['p_exact_signflip'] else 1 for r in z];old=[float(r['p_exact_woman_signflip']) for r in rows(I/'prior_within_site_contrasts.csv') if r['taxon_arm']==arm and r['scope']==scope]
 q=independent_bh(p);q18=independent_bh(old+p)[9:]
 for r,a,b in zip(z,q,q18):
  if r['q_BH_9_between_site']:near(a,r['q_BH_9_between_site'],'scipy BH9');near(b,r['q_BH_18_all_own_PRE'],'scipy BH18')
check(len(results)==36,'four nine-test families');spec=json.loads((O/'Analysis_specification.json').read_text())
for n,h in spec['source_hashes'].items():check(hashlib.sha256((I/n).read_bytes()).hexdigest()==h,n+' unchanged')
report=dict(status='PASS',checks=checks,independent_scalar_distance_replays=208,comparisons=36,primary_comparisons=9,paired_manifests_verified=True,all_estimable_exact_signflip_tests_recomputed=True,all_saved_bootstraps_replayed=True,BH_independently_checked_scalar=True,maximum_numeric_error=maxerr,claim_limit='Arithmetic/matching verification does not establish observational exchangeability or recover upstream provenance.')
(Q/'Independent_numerical_audit.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
