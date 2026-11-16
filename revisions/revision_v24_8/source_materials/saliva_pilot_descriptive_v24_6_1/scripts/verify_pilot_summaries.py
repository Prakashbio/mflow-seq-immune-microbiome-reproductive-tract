"""Independent read-only FCS decoding and verification of supplied pilot summaries.
No gates, compensation, donor pairing or antibody-specificity inference is applied.
"""
from pathlib import Path
from collections import Counter
import csv,json,hashlib,re,struct,argparse,datetime
import numpy as np
R=Path(__file__).resolve().parents[1];p=argparse.ArgumentParser();p.add_argument('--source',type=Path,default=R/'Prepared_repository_update_v24_6_2026-11-12/facs_control_pilot_display_v23_2/source_fcs');p.add_argument('--out',type=Path,default=R/'06_checks/V24_6_final_comment_reaudit_2026-11-12/saliva_pilot_independent');p.add_argument('--provided-folder',type=Path);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
def fcs(path):
 raw=path.read_bytes();assert raw[:3]==b'FCS'
 start,end=[int(raw[i:i+8]) for i in [10,18]];txt=raw[start:end+1].decode('latin-1');delimiter=txt[0];tokens=txt[1:].split(delimiter)
 if tokens[-1]=='':tokens.pop()
 assert len(tokens)%2==0
 m=dict(zip(tokens[::2],(x.strip() for x in tokens[1::2])));n=int(m['$PAR']);total=int(m['$TOT']);assert m['$DATATYPE']=='F';assert all(m[f'$P{i}B']=='32' for i in range(1,n+1));offset=int(m['$BEGINDATA']);size=n*total*4;assert int(m['$ENDDATA'])-offset+1==size
 endian='>' if m['$BYTEORD']=='4,3,2,1' else '<' if m['$BYTEORD']=='1,2,3,4' else None;assert endian
 # Independent struct check of selected records against the bulk NumPy view.
 data=np.frombuffer(raw[offset:offset+size],dtype=endian+'f4').reshape(total,n).astype('float64');assert np.isfinite(data).all()
 for row in [0,total//2,total-1]:assert np.array_equal(data[row],struct.unpack_from(endian+'f'*n,raw,offset+row*n*4))
 return raw,m,data,[m[f'$P{i}N'] for i in range(1,n+1)]
def outcsv(name,rows):
 with (a.out/name).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
summary=[];settings=[];files={};checks=[]
for path in sorted(a.source.rglob('*.fcs')):
 raw,m,data,names=fcs(path);digest=hashlib.sha256(raw).hexdigest();files[digest]=(path,m,data,names);spill=m.get('$SPILL',m.get('SPILL',''));identity=None
 if spill:
  vals=spill.split(',');k=int(vals[0]);mat=np.array(list(map(float,vals[k+1:]))).reshape(k,k);identity=bool(np.array_equal(mat,np.eye(k)))
 base=dict(date=m['$DATE'],tube=m.get('TUBE NAME',''),FCS_SHA256=digest,events=len(data),start_time=m.get('$BTIM',''),end_time=m.get('$ETIM',''),threshold=m.get('THRESHOLD',''),SPILL_identity=identity,APPLY_COMP=m.get('APPLY COMPENSATION','not-recorded'),source=str(path.relative_to(a.source)))
 for i,name in enumerate(names,1):
  settings.append(dict(**base,channel=name,voltage=m.get(f'$P{i}V',''),gain=m.get(f'$P{i}G',''),range=m.get(f'$P{i}R',''),amplification=m.get(f'$P{i}E','')))
  if name not in {'2-488-530/30-A','2-488-616/23-A','2-488-585/42-A','1-633-660/20-A','3-406-450/40-A'}:continue
  x=data[:,i-1];qq=np.quantile(x,[.05,.25,.5,.75,.95]);summary.append(dict(**base,channel=name,median=float(qq[2]),q05=float(qq[0]),q25=float(qq[1]),q75=float(qq[3]),q95=float(qq[4]),fraction_nonpositive=float(np.mean(x<=0)),max=float(max(x))))
outcsv('Independently_recomputed_channel_summary.csv',summary);outcsv('All_acquisition_detector_settings.csv',settings)
provided=list(csv.DictReader((a.provided_folder/'pilot_fcs_channel_summary.csv').open())) if a.provided_folder else [];lookup={(x['FCS_SHA256'],x['channel']):x for x in summary};maxerr=0.;mismatches=[];precision=[]
for row in provided:
 key=(row['FCS_SHA256'],row['channel']);assert key in lookup;z=lookup[key]
 for col in ['median','q05','q25','q75','q95','fraction_nonpositive','max','events']:
  err=abs(float(row[col])-float(z[col]));maxerr=max(maxerr,err)
  if err>1e-8:
   path,meta,data,names=files[row['FCS_SHA256']]
   if col in ['q05','q25','q75','q95']:
    qv={'q05':.05,'q25':.25,'q75':.75,'q95':.95}[col];reproduced=float(np.quantile(data[:,names.index(row['channel'])].astype('float32'),np.float32(qv)))
    if abs(reproduced-float(row[col]))<1e-7:precision.append(dict(FCS_SHA256=row['FCS_SHA256'],channel=row['channel'],field=col,supplied_float32_quantile=row[col],independent_float64_quantile=z[col],absolute_difference=err));continue
   mismatches.append(dict(key=str(key),field=col,provided=row[col],computed=z[col],absolute_error=err))
 for col in ['date','tube','start_time','end_time','threshold','SPILL_identity','APPLY_COMP']:
  if str(row[col])!=str(z[col]):mismatches.append(dict(key=str(key),field=col,provided=row[col],computed=z[col],absolute_error='text'))
assert len(lookup)==167 and (not a.provided_folder or len(provided)==167);assert not mismatches,mismatches
if precision:outcsv('Float_precision_differences.csv',precision)
# Check intended detector settings across each matched-by-index comparison before interpreting medians.
matched=[];timepairs=[]
def one(date,tube):
 hits=[v for v in files.values() if v[1]['$DATE']==date and v[1].get('TUBE NAME')==tube];assert len(hits)==1;return hits[0]
def setting(m,names,channel):
 j=names.index(channel)+1;return {k:m.get(k,'') for k in [f'$P{j}V',f'$P{j}G',f'$P{j}R',f'$P{j}E','THRESHOLD','$SPILL','APPLY COMPENSATION']}
for i in range(1,5):
 u=one('13-JUN-2023',f'Saliva unstained {i}-4')
 for label,channel,suffix in [('IgA','2-488-585/42-A','2'),('IgM','1-633-660/20-A','3')]:
  l=one('13-JUN-2023',f'Saliva {label}+'+('PE' if label=='IgA' else 'APC')+f' {i}-{suffix}');lm=float(np.median(l[2][:,l[3].index(channel)]));um=float(np.median(u[2][:,u[3].index(channel)]));ls=setting(l[1],l[3],channel);us=setting(u[1],u[3],channel)
  matched.append(dict(acquisition_index=i,label=label,channel=channel,labelled_median=lm,unstained_median=um,difference=lm-um,labelled_events=len(l[2]),unstained_events=len(u[2]),detector_and_trigger_settings_match=ls==us,labelled_settings=json.dumps(ls),unstained_settings=json.dumps(us),labelled_FCS_SHA256=hashlib.sha256(l[0].read_bytes()).hexdigest(),unstained_FCS_SHA256=hashlib.sha256(u[0].read_bytes()).hexdigest()))
 x=one('20-JUN-2023',f'Saliva Syto9+PI {i}-1');y=one('20-JUN-2023',f'Saliva IgA+PE {i}-2');t=lambda z:datetime.datetime.strptime(z[1]['$BTIM'],'%H:%M:%S');gap=(t(y)-t(x)).total_seconds()/60;ch='3-406-450/40-A'
 timepairs.append(dict(acquisition_index=i,SYTO9_PI_start=x[1]['$BTIM'],DAPI_IgA_start=y[1]['$BTIM'],interval_between_starts_min=gap,SYTO9_PI_events=len(x[2]),DAPI_IgA_events=len(y[2]),violet_SYTO9_PI_median=float(np.median(x[2][:,x[3].index(ch)])),violet_DAPI_IgA_median=float(np.median(y[2][:,y[3].index(ch)])),detector_and_trigger_settings_match=setting(x[1],x[3],ch)==setting(y[1],y[3],ch)))
outcsv('Independent_June13_channel_contrasts.csv',matched);outcsv('Independent_June20_candidate_pairs.csv',timepairs)
for row in (csv.DictReader((a.provided_folder/'2023_06_13_saliva_antibody_channel_contrasts.csv').open()) if a.provided_folder else []):
 z=next(x for x in matched if str(x['acquisition_index'])==row['acquisition_index'] and x['label']==row['label']);assert abs(round(z['difference'],2)-float(row['median_difference_all_events']))<1e-8
for row in (csv.DictReader((a.provided_folder/'2023_06_20_saliva_candidate_pairs.csv').open()) if a.provided_folder else []):
 z=next(x for x in timepairs if str(x['acquisition_index'])==row['acquisition_index']);assert abs(round(z['interval_between_starts_min'],2)-float(row['interval_between_starts_min']))<1e-8
report=dict(status='PASS',FCS_files=len(files),fluorescence_channel_rows=len(summary),supplied_summary_numeric_cells_verified=len(provided)*8,max_absolute_difference=maxerr,float32_percentile_differences_reproduced=len(precision),all_medians_and_counts_match=True if a.provided_folder else None,antibody_contrasts=len(matched),all_June13_compared_settings_identical=all(x['detector_and_trigger_settings_match'] for x in matched),all_June20_compared_settings_identical=all(x['detector_and_trigger_settings_match'] for x in timepairs),acquisitions_by_date=dict(Counter(v[1]['$DATE'] for v in files.values())),additional_gating_or_compensation_applied=False,scope='All recorded event-array summaries. Numbered acquisition labels are not verified donor pairing; medians do not validate antibody specificity or DAPI agreement.')
(a.out/'Independent_pilot_verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2));print('SETTINGS CONTRASTS',[(x['label'],x['acquisition_index'],x['detector_and_trigger_settings_match']) for x in matched]);print('JUNE20',timepairs)
