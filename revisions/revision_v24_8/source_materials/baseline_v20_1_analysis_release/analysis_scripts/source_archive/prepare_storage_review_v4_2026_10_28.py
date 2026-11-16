"""Reviewer-requested descriptive storage review; no preservation/equivalence test."""
from pathlib import Path
import os
import argparse,hashlib,json,shutil
import numpy as np
import pandas as pd
R=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path);p.add_argument('--output',type=Path,default=R/'storage_integrated_working_v4_2026-10-28');a=p.parse_args();O=a.output
I=O/'storage_inputs';D=O/'storage_statistics';I.mkdir(parents=True,exist_ok=True);D.mkdir(exist_ok=True)
sha=lambda q:hashlib.sha256(q.read_bytes()).hexdigest()
if a.inputs:
 for q in a.inputs.iterdir():
  if q.is_file() and q.resolve()!=(I/q.name).resolve():shutil.copy2(q,I/q.name)
else:
 T=R/'03_analysis/Taxon_scores_PRE_review_2026-10-28/data'
 for n in ['baseline_39_source_counts.csv','sensitivity_41_source_counts.csv','source_metadata_original_codes.csv','original_taxonomy.csv']:shutil.copy2(T/n,I/n)
 manifest=R/'02_evidence/Source_records_2026-10-20/Source_record_manifest_REVIEW_ONLY.csv';q=pd.read_csv(manifest)
 q=q[['SampleID','Participant','Location','Antibody','Storage','workbook_row','workbook_patient_cell','workbook_site_cell','profile_sheet','profile_header_cell','existing_processed_heading','report_sha256']].rename(columns={'SampleID':'SourceCode'})
 q.to_csv(I/'Storage_workbook_source_links.csv',index=False)
 import openpyxl
 native=Path(os.environ.get('MFLOW_NATIVE_WORKBOOK', 'Prakash FACS study final version of microbiome sequencing. Improved pipline 07.04.25.xlsx'));before=sha(native)
 wb=openpyxl.load_workbook(native,data_only=True,read_only=True);sheet=wb['General table'];snapshot=[]
 for row in q.itertuples():
  header='A1' if int(row.workbook_row)<29 else 'A29';recorded=str(sheet[header].value).split()[0].upper();patient=str(sheet[row.workbook_patient_cell.split('!')[-1]].value);site=str(sheet[row.workbook_site_cell.split('!')[-1]].value)
  assert recorded==row.Storage and int(patient.removeprefix('HUT'))==row.Participant
  assert {'Vaginal Swab':'VAGINA','Cervix':'CERVIX','Endometrium':'ENDOMETRIUM'}[site]==row.Location
  snapshot.append(dict(SourceCode=row.SourceCode,Storage=row.Storage,workbook_row=row.workbook_row,storage_header_cell=header,storage_header_value=str(sheet[header].value),patient_cell=row.workbook_patient_cell,patient_value=patient,site_cell=row.workbook_site_cell,site_value=site,passed=True))
 wb.close();assert sha(native)==before
 pd.DataFrame(snapshot).to_csv(I/'Storage_native_workbook_cell_verification.csv',index=False)
 shutil.copy2(R/'02_evidence/Original_clinical_methods_2026-10-26/Merli_original_methods_extracted.txt',I/'Merli_original_methods_extracted.txt')
 shutil.copy2(R/'02_evidence/Native_storage_trace_2026-10-28/HUT099_native_storage_code_trace.csv',I/'HUT099_native_storage_code_trace.csv')
 (I/'Input_origins.json').write_text(json.dumps(dict(source_manifest=str(manifest),source_manifest_sha256=sha(manifest),native_workbook=str(native),native_workbook_sha256=before,native_storage_patient_site_checks=149,storage_fields='Values independently rechecked against General table storage section, patient and site cells; not inferred from suffix or derived table headings.',originals_unchanged=True),indent=2)+'\n')
hashes={q.name:sha(q) for q in I.iterdir() if q.is_file()}
m=pd.read_csv(I/'source_metadata_original_codes.csv');links=pd.read_csv(I/'Storage_workbook_source_links.csv')
assert len(m)==len(links)==149 and links.SourceCode.is_unique and m.SourceCode.is_unique
m=m.merge(links[['SourceCode','Participant','Location','Storage','workbook_row']],on=['SourceCode','Participant','Location'],validate='one_to_one')
assert len(m)==149 and m.Storage.isin(['FRESH','FROZEN']).all() and m.Participant.nunique()==15 and m.SpecimenKey.nunique()==40
assert m.groupby('SpecimenKey').Storage.nunique().max()==1
assert list(m[(m.Participant==99)&(m.Location=='ENDOMETRIUM')].groupby('Storage').size())==[4,4]
repeats=m.groupby(['Participant','Location']).Storage.nunique();assert repeats[repeats>1].index.tolist()==[(99,'ENDOMETRIUM')]
md=m.copy();md['quantity_note']='One retained sequencing library; a preparation/fraction observation, not an independent woman or FACS percentage.';md.to_csv(D/'Storage_library_manifest.csv',index=False)
tax=pd.read_csv(I/'original_taxonomy.csv').set_index('Taxon');metrics=[];relrows=[];means=[];summaries=[];pairs=[];pairsmetrics=[]
for arm in ['baseline_39','sensitivity_41']:
 x=pd.read_csv(I/f'{arm}_source_counts.csv',index_col=0);assert x.columns.tolist()==m.SourceCode.tolist() and (x.sum()>0).all()
 ra=x.div(x.sum(axis=0),axis=1);lacto=tax.loc[x.index,'Genus'].eq('Lactobacillus');lactotaxa=x.index[lacto].tolist();assert len(lactotaxa)==5
 for row in m.itertuples():
  v=ra[row.SourceCode];pos=v[v>0];H=float(-(pos*np.log(pos)).sum());n=len(pos)
  metrics.append(dict(taxon_arm=arm,SourceCode=row.SourceCode,Participant=row.Participant,SpecimenKey=row.SpecimenKey,Location=row.Location,Fraction=row.Fraction,Storage=row.Storage,Shannon=H,Simpson_1_minus_D=float(1-(v*v).sum()),Richness=n,Pielou_evenness=H/np.log(n) if n>1 else np.nan,Lactobacillus_genus_pct=float(100*v[lacto].sum()),retained_count_sum=int(x[row.SourceCode].sum())))
  for name,val in v.items():relrows.append(dict(taxon_arm=arm,SourceCode=row.SourceCode,Participant=row.Participant,SpecimenKey=row.SpecimenKey,Location=row.Location,Fraction=row.Fraction,Storage=row.Storage,Taxon=name,relative_abundance=float(val)))
 # Averaging makes the contribution rule explicit if future source records
 # contain more than one preparation in a woman/site/condition/fraction cell.
 lib=pd.DataFrame([z for z in metrics if z['taxon_arm']==arm]);woman=lib.groupby(['Participant','Location','Fraction','Storage'])
 for (site,frac,storage),g in lib.groupby(['Location','Fraction','Storage']):
  wm=g.groupby('Participant')[['Shannon','Simpson_1_minus_D','Richness','Pielou_evenness','Lactobacillus_genus_pct']].mean()
  for metric in wm:
   v=wm[metric].dropna();summaries.append(dict(taxon_arm=arm,Location=site,Fraction=frac,Storage=storage,metric=metric,n_libraries=len(g),n_preparations=g.SpecimenKey.nunique(),n_women=len(wm),n_nonmissing=len(v),n_missing=len(wm)-len(v),mean=float(v.mean()),SD=float(v.std(ddof=1)),median=float(v.median()),Q1=float(v.quantile(.25)),Q3=float(v.quantile(.75)),minimum=float(v.min()),maximum=float(v.max()),unit='Equal-woman descriptive summary within site/fraction/recorded storage; no causal or equivalence inference.'))
  wra=[]
  for _,person in g.groupby('Participant'):wra.append(ra[person.SourceCode].mean(axis=1))
  avg=pd.concat(wra,axis=1).mean(axis=1)
  for t,v in avg.items():means.append(dict(taxon_arm=arm,Location=site,Fraction=frac,Storage=storage,Taxon=t,mean_relative_abundance=float(v),n_women=len(wra),n_libraries=len(g)))
 z=m[(m.Participant==99)&(m.Location=='ENDOMETRIUM')]
 for frac,g in z.groupby('Fraction'):
  f=g.loc[g.Storage=='FRESH','SourceCode'].item();c=g.loc[g.Storage=='FROZEN','SourceCode'].item();fv=ra[f];cv=ra[c]
  pairs.append(dict(taxon_arm=arm,Participant=99,Location='ENDOMETRIUM',Fraction=frac,fresh_SourceCode=f,frozen_SourceCode=c,Bray_Curtis_relative_abundance=float(np.abs(fv-cv).sum()/(fv.sum()+cv.sum())),n_women=1,n_recorded_preparations=2,interpretation='Same woman and two recorded preparations; not verified split aliquots, replicated storage effect or an equivalence test.'))
  for metric in ['Shannon','Lactobacillus_genus_pct']:
   lf=lib.loc[lib.SourceCode==f,metric].item();lc=lib.loc[lib.SourceCode==c,metric].item();pairsmetrics.append(dict(taxon_arm=arm,Fraction=frac,metric=metric,fresh_value=lf,frozen_value=lc,frozen_minus_fresh=lc-lf,n_women=1,interpretation='Descriptive values only; fractions are not independent storage replications.'))
for n,rows in [('Storage_library_metrics',metrics),('Storage_library_taxon_relative',relrows),('Storage_equal_woman_group_means',means),('Storage_equal_woman_descriptive_summaries',summaries),('Participant99_two_preparation_distances',pairs),('Participant99_two_preparation_metric_values',pairsmetrics)]:pd.DataFrame(rows).to_csv(D/f'{n}.csv',index=False)
den=m.groupby(['Location','Fraction','Storage'],as_index=False).agg(n_libraries=('SourceCode','size'),n_preparations=('SpecimenKey','nunique'),n_women=('Participant','nunique'));den.to_csv(D/'Storage_group_denominators.csv',index=False)
prep=m.drop_duplicates('SpecimenKey');prep.groupby(['Location','Storage'],as_index=False).agg(n_preparations=('SpecimenKey','size'),n_women=('Participant','nunique')).to_csv(D/'Storage_preparation_denominators.csv',index=False)
overlap=[]
for (site,frac),g in m.groupby(['Location','Fraction']):
 fresh=set(g.loc[g.Storage=='FRESH','Participant']);frozen=set(g.loc[g.Storage=='FROZEN','Participant']);overlap.append(dict(Location=site,Fraction=frac,n_fresh_women=len(fresh),n_frozen_women=len(frozen),n_union_women=len(fresh|frozen),n_women_in_both_conditions=len(fresh&frozen),participants_in_both=';'.join(map(str,sorted(fresh&frozen)))))
pd.DataFrame(overlap).to_csv(D/'Storage_between_condition_woman_overlap.csv',index=False)
assert all(sha(I/n)==h for n,h in hashes.items())
result=dict(status='PASS_DESCRIPTIVE_PREPARATION',input_hashes=hashes,libraries=149,women=15,preparations=40,group_rows=len(den),frozen_cervical_PRE_women=2,source_codes_frozen_cervical_PRE=m.loc[(m.Location=='CERVIX')&(m.Fraction=='PRE')&(m.Storage=='FROZEN'),'SourceCode'].tolist(),only_cross_condition_woman_site='Participant99,endometrium',randomized_storage_design_documented=False,verified_split_aliquot_pairs=False,equivalence_test_performed=False,hypothesis_tests_performed=False,FACS_percentage_link_established=False,limits='Conditional on named-source/workbook mapping. Historical derived-heading correspondence and upstream identity require author adjudication. Storage is mostly between women; variation cannot be attributed specifically to freezing. Cohort composition, menstrual phase, processing and storage history are not controlled by these descriptions.')
(O/'Storage_analysis_provenance.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='input_hashes'},indent=2))
