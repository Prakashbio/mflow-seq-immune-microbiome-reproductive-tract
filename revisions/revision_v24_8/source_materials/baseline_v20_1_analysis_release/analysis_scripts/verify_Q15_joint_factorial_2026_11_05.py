"""Independent NumPy/pandas verification of the supplemental joint analysis.

Checks numerical implementation, not the truth of residual exchangeability.
Usage: python verify_Q15_joint_factorial_2026_11_05.py ANALYSIS_DIRECTORY [--report-out REPORT_JSON]
The analysis directory is read-only. The optional report must be outside it.
"""
import sys, json, hashlib, argparse
from pathlib import Path
import numpy as np
import pandas as pd

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('analysis_directory',type=Path)
parser.add_argument('--report-out',type=Path,help='Optional JSON report outside the analysis directory; default prints to stdout only')
args=parser.parse_args()
o=args.analysis_directory.resolve()
report_out=args.report_out.resolve() if args.report_out else None
if report_out is not None and (report_out==o or o in report_out.parents):
    parser.error('--report-out must be outside the read-only analysis directory')
if report_out is not None and report_out.exists():
    parser.error('--report-out already exists; choose a new filename to preserve prior reports')
res=o/'results'; diag=o/'diagnostics'
tab=pd.read_csv(res/'Joint_factorial_marginal_tests_and_intervals.csv')
m=pd.read_csv(o/'inputs/Metadata_original_codes.csv',dtype={'Participant':str,'SourceSuffix':str}).fillna('')
P=pd.read_csv(diag/'Hierarchical_permutation_indices_1based.csv').to_numpy(int)-1
null=pd.read_csv(diag/'All_marginal_null_F_statistics.csv')
boots=pd.read_csv(diag/'Woman_bootstrap_indices_1based.csv').to_numpy(int)-1
key=pd.read_csv(diag/'Bootstrap_woman_key.csv',dtype={'Participant':str}).Participant.to_numpy()
checks=[]
def check(label, ok, detail=None):
    checks.append({'check':label,'passed':bool(ok),'detail':detail})
    if not ok: raise AssertionError((label,detail))
def hat(x):
    u,s,v=np.linalg.svd(x,full_matrices=False)
    rank=int((s>max(x.shape)*np.finfo(float).eps*s[0]).sum())
    return u[:,:rank]@u[:,:rank].T,rank
def correction(p,method):
    p=np.asarray(p); ix=np.argsort(p); out=np.empty_like(p); n=len(p)
    out[ix]=np.minimum(1,np.maximum.accumulate(p[ix]*(n-np.arange(n)))) if method=='holm' else np.minimum(1,np.minimum.accumulate((p[ix]*n/(np.arange(n)+1))[::-1])[::-1])
    return out
max_null_error=0.
for tag,rows in tab.groupby('analysis',sort=False):
    zz=pd.read_csv(res/(tag+'_profile_manifest.csv'),dtype={'Participant':str})
    yfile=pd.read_csv(res/(tag+'_relative_abundance_profiles.csv'))
    y=yfile.iloc[:,1:].to_numpy(float)
    taxa=int(rows.taxon_membership.iloc[0]); counts=pd.read_csv(o/'inputs'/f'Counts_{taxa}_original_codes.csv',index_col=0)
    relative=counts/counts.sum(axis=0)
    rebuilt=np.stack([relative[s.split(';')].mean(axis=1).to_numpy() for s in zz.source_libraries])
    check(tag+' normalized profiles independently reconstructed',np.allclose(y,rebuilt,atol=1e-14,rtol=0),float(np.max(abs(y-rebuilt))))
    used=[s for v in zz.source_libraries for s in v.split(';')]
    check(tag+' source library count',len(used)==int(rows.underlying_libraries.iloc[0]) and len(set(used))==len(used),len(used))
    if tag.endswith('mean_both'):
        double=zz[zz.source_preparations==2]
        check(tag+' HUT099 both preparations equally represented',len(double)==4 and set(double.Participant)=={'99'} and set(double.Location)=={'ENDOMETRIUM'})
    D=np.abs(y[:,None,:]-y[None,:,:]).sum(axis=2)/(y[:,None,:]+y[None,:,:]).sum(axis=2)
    savedD=pd.read_csv(diag/(tag+'_BrayCurtis.csv')).to_numpy(float)
    check(tag+' independently calculated Bray Curtis',np.allclose(D,savedD,atol=1e-12,rtol=0),float(np.max(abs(D-savedD))))
    n=len(zz); J=np.eye(n)-np.ones((n,n))/n; G=-.5*J@(D*D)@J
    full=pd.read_csv(diag/(tag+'_full_design.csv')).to_numpy(float); Hf,rank=hat(full); Rf=np.eye(n)-Hf
    check(tag+' common full rank and residual df',rank==25 and n-rank==95)
    # All permutations preserve independent women and both covariance relations.
    if tag==tab.analysis.iloc[0]:
        check('every row is a permutation',np.all(np.sort(P,axis=1)==np.arange(n)))
        person=zz.Participant.to_numpy(); unit=zz.AnatomicalUnit.to_numpy(); wf=(zz.Participant+'_'+zz.Fraction).to_numpy()
        check('no permutation crosses a woman',np.all(person[P]==person[None,:]))
        for label,v in [('anatomical unit',unit),('woman fraction',wf)]:
            rel=v[:,None]==v[None,:]
            check('all permutations preserve '+label,all(np.array_equal(rel[np.ix_(p,p)],rel) for p in P),len(P))
    for r in rows.itertuples():
        red=pd.read_csv(diag/(tag+'_'+r.term+'_reduced_design.csv')).to_numpy(float)
        Hr,rr=hat(red); Ht=Hf-Hr
        ss=np.sum(G*Ht); rss=np.sum(G*Rf); f=(ss/(rank-rr))/(rss/(n-rank))
        check(tag+' '+r.term+' SS F and R2',np.allclose([ss,f,ss/np.trace(G)],[r.marginal_SS,r.pseudo_F,r.term_R2_total],atol=1e-9,rtol=1e-9))
        nv=null[tag+'__'+r.term].to_numpy()
        p=(1+(nv>=r.pseudo_F-np.sqrt(np.finfo(float).eps)).sum())/(len(P)+1)
        check(tag+' '+r.term+' permutation p',abs(p-r.p_raw)<1e-14)
        # Independently reconstruct selected null statistics from the residual
        # distance geometry. Neither vegan nor its fitted object is used here.
        Rr=np.eye(n)-Hr; Gres=Rr@G@Rr
        indices=np.unique(np.linspace(0,len(P)-1,min(199,len(P))).astype(int))
        calc=[]
        for k in indices:
            gp=Gres[np.ix_(P[k],P[k])]
            calc.append(np.sum(gp*Ht)/(rank-rr)/(np.sum(gp*Rf)/(n-rank)))
        err=float(np.max(abs(np.array(calc)-nv[indices]))); max_null_error=max(max_null_error,err)
        check(tag+' '+r.term+' independent null statistics',err<1e-8,{'null_permutations_recomputed':len(indices),'max_error':err})
    check(tag+' Holm independently reproduced',np.allclose(correction(rows.p_raw,'holm'),rows.p_Holm_three_terms,atol=1e-14))
    check(tag+' BH independently reproduced',np.allclose(correction(rows.p_raw,'bh'),rows.p_BH_three_terms,atol=1e-14))
    bv=pd.read_csv(res/(tag+'_R2_woman_bootstrap_draws.csv'))
    for r in rows.itertuples():
        v=bv[r.term].dropna().to_numpy(); ci=np.quantile(v,[.025,.975])
        check(tag+' '+r.term+' bootstrap percentile interval',len(v)==r.valid_bootstrap_draws and np.allclose(ci,[r.R2_CI95_low,r.R2_CI95_high],atol=1e-12))
    # Rebuild 50 cluster draws using an independently coded sum-contrast design.
    sites=['VAGINA','CERVIX','ENDOMETRIUM']; fracs=['PRE','IgA','IgM','IgG']
    def sumcodes(values,levels):
        code=np.vstack([np.eye(len(levels)-1),-np.ones((1,len(levels)-1))])
        return code[[levels.index(v) for v in values]]
    for b in range(min(50,len(boots))):
        ix=np.concatenate([np.flatnonzero(zz.Participant.to_numpy()==key[k]) for k in boots[b]])
        sizes=[int((zz.Participant==key[k]).sum()) for k in boots[b]]
        Z=np.repeat(np.eye(len(key)),sizes,axis=0)
        loc=sumcodes(zz.Location.iloc[ix].tolist(),sites); fra=sumcodes(zz.Fraction.iloc[ix].tolist(),fracs)
        inter=(fra[:,:,None]*loc[:,None,:]).reshape(len(ix),-1)
        terms=[loc,fra,inter]; ff=np.column_stack([Z,*terms]); hh,rrank=hat(ff)
        jj=np.eye(len(ix))-np.ones((len(ix),len(ix)))/len(ix); gg=-.5*jj@(D[np.ix_(ix,ix)]**2)@jj
        if rrank!=25: continue
        vals=[]
        for t in range(3):
            hr,_=hat(np.column_stack([Z,*[v for k,v in enumerate(terms) if k!=t]]))
            vals.append(np.sum(gg*(hh-hr))/np.trace(gg))
        expected=bv.loc[b,['Site','Fraction','Site_by_fraction']].to_numpy(float)
        check(tag+' independent cluster draw '+str(b+1),np.allclose(vals,expected,atol=1e-9,rtol=1e-9))
disp=pd.read_csv(res/'Table_S2R_10_joint_subset_dispersion.csv')
for tag,g in disp.groupby('analysis',sort=False):
    check(tag+' dispersion Holm adjustment',np.allclose(correction(g.p_raw,'holm'),g.p_Holm_two_dispersion_tests,atol=1e-14))
    for r in g.itertuples():
        dd=pd.read_csv(diag/(tag+'_'+r.grouping+'_dispersion_distances.csv'))
        values=dd.distance_to_spatial_median.to_numpy()
        groups=dd.group.to_numpy(); levels=np.unique(groups)
        mu=np.array([values[groups==v].mean() for v in levels])
        ng=np.array([(groups==v).sum() for v in levels])
        means=np.array([mu[list(levels).index(v)] for v in groups])
        residual=values-means
        f=(np.sum(ng*(mu-values.mean())**2)/(len(levels)-1))/(np.sum(residual**2)/(len(values)-len(levels)))
        yp=residual[P]
        gm=yp.mean(axis=1)
        mm=np.stack([yp[:,groups==v].mean(axis=1) for v in levels],axis=1)
        mss=np.sum(ng[None,:]*(mm-gm[:,None])**2,axis=1)
        rss=np.sum(yp**2,axis=1)-np.sum(ng[None,:]*mm**2,axis=1)
        ff=(mss/(len(levels)-1))/(rss/(len(values)-len(levels)))
        pv=(1+(ff>=f-np.sqrt(np.finfo(float).eps)).sum())/(len(P)+1)
        check(tag+' '+r.grouping+' dispersion F and all permutation p',np.allclose([f,pv],[r.F,r.p_raw],rtol=1e-10,atol=1e-12))
result={'status':'PASS','checks_passed':len(checks),'permutations':len(P),'independent_null_statistics_recomputed':12*min(199,len(P)),'maximum_null_F_error':max_null_error,'checks':checks,'scope':'Numerical verification only. Does not establish biological validity, residual exchangeability, or equal multivariate dispersion.'}
if report_out is not None:
    report_out.parent.mkdir(parents=True,exist_ok=True)
    with report_out.open('x') as report:
        report.write(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='checks'},indent=2))
