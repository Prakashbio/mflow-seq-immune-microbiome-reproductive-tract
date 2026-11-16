"""Source-filter audit. No new contamination classification or source inference."""
from pathlib import Path
import argparse,ast,hashlib,json,shutil
import numpy as np
import pandas as pd
from openpyxl import load_workbook

R=Path(__file__).resolve().parents[1];p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path);p.add_argument('--output',type=Path);a=p.parse_args()
O=a.output or R/'source_filter_integrated_working_v5_2026-10-28';I=O/'source_filter_inputs';D=O/'source_filter_statistics';I.mkdir(parents=True,exist_ok=True);D.mkdir(exist_ok=True)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
if a.inputs:
 for f in a.inputs.iterdir():
  if f.is_file() and f.resolve()!=(I/f.name).resolve():shutil.copy2(f,I/f.name)
else:
 root=R.parent;v4=R/'storage_integrated_working_v4_2026-10-28';e=R/'02_evidence/Source_records_2026-10-20';g=R/'02_evidence/Public_repository_access_2026-10-26/GitHub_snapshot'
 sources={
  'Counts_268_original_codes.csv':root/'mFLOW-seq_Menuscript/Zenodo home/ReadsTable_P.csv',
  'Metadata_original_codes.csv':v4/'storage_inputs/source_metadata_original_codes.csv',
  'Counts_41_original_codes.csv':v4/'storage_inputs/sensitivity_41_source_counts.csv',
  'Counts_39_original_codes.csv':v4/'storage_inputs/baseline_39_source_counts.csv',
  'Workbook_link_manifest.csv':e/'Source_record_manifest_REVIEW_ONLY.csv',
  'Recovered_reference_178.xlsx':root/'mFLOW-seq_Menuscript/Menuscript Final/Addtional data files/ContaminantReference_178taxa_MerliSaare_2026.xlsx',
  'Original_sfig16_script.py':g/'scripts/sfig16_contamination_analysis.py',
  'Original_seq_clean.csv':g/'data/seq_clean.csv',
  'Original_otu_clean.csv':g/'data/otu_clean.csv',
  'Original_otu_clean_decontaminated.csv':g/'data/otu_clean_decontaminated.csv'}
 for name,src in sources.items():shutil.copy2(src,I/name)
 (I/'Input_origins.json').write_text(json.dumps({n:dict(original_path=str(f),sha256=sha(f)) for n,f in sources.items()},indent=2)+'\n')
x=pd.read_csv(I/'Counts_268_original_codes.csv',index_col=0);x41=pd.read_csv(I/'Counts_41_original_codes.csv',index_col=0);x39=pd.read_csv(I/'Counts_39_original_codes.csv',index_col=0)
m=pd.read_csv(I/'Metadata_original_codes.csv');links=pd.read_csv(I/'Workbook_link_manifest.csv').set_index('SampleID')
assert x.shape==(268,149) and list(x.columns)==m.SourceCode.tolist()==list(x41.columns)==list(x39.columns)
assert np.array_equal(x.loc[x41.index],x41) and np.array_equal(x41.loc[x39.index],x39)
assert set(x41.index)-set(x39.index)=={'Burkholderia cepacia','Streptococcus pneumoniae'}
groups=m.groupby(['Location','Fraction'],sort=True)
ra=x.div(x.sum(),axis=1);ra_workbook=x.div(links.loc[x.columns,'workbook_reads_numeric'].to_numpy(),axis=1)
rules=[('0.5% in 2 libraries',.005,2),('0.5% in 3 libraries',.005,3),('0.5% in 4 libraries',.005,4),('10% in 2 libraries',.10,2)]
thresholds=[];passes=pd.DataFrame(index=x.index)
for (site,frac),md in groups:
 ids=md.SourceCode.tolist()
 for label,threshold,nmin in rules:
  count=(ra[ids]>=threshold).sum(axis=1);kept=count>=nmin
  thresholds.append(dict(Location=site,Fraction=frac,rule=label,minimum_relative_abundance=threshold,minimum_libraries=nmin,n_group_libraries=len(ids),n_taxa_passing=int(kept.sum())))
  if nmin==2 and threshold==.005:passes[site+'_'+frac]=kept
selected=passes.any(axis=1);selected_workbook=np.logical_or.reduce([((ra_workbook[md.SourceCode]>=.005).sum(axis=1)>=2).to_numpy() for _,md in groups])
assert set(x.index[selected])==set(x41.index)==set(x.index[selected_workbook])
pd.DataFrame(thresholds).to_csv(D/'Table_S2C_3_threshold_counts.csv',index=False)
inventory=pd.DataFrame(dict(Taxon=x.index,in_41=x.index.isin(x41.index),in_39=x.index.isin(x39.index),total_exported_counts=x.sum(axis=1),libraries_detected=(x>0).sum(axis=1),passes_reconstructed_rule=selected,groups_passing=passes.sum(axis=1)))
inventory['record_status']=np.where(inventory.in_39,'Retained submitted baseline',np.where(inventory.in_41,'Historical manual exclusion; contaminant status not established','Below reconstructed analytical selection rule; contaminant status not established'))
inventory.to_csv(D/'Table_S2C_2_exported_taxon_membership.csv',index=False)
q=m[['SourceCode','Participant','Location','Fraction','SpecimenKey']].copy()
q['exported_268_counts']=x.sum().to_numpy();q['retained_41_counts']=x41.sum().to_numpy();q['retained_39_counts']=x39.sum().to_numpy()
q['discarded_268_to_41_pct']=100*(1-q.retained_41_counts/q.exported_268_counts)
q['discarded_41_to_39_pct']=100*(1-q.retained_39_counts/q.retained_41_counts)
q.to_csv(D/'Filtering_library_disposition.csv',index=False)
w=load_workbook(I/'Recovered_reference_178.xlsx',read_only=True,data_only=True);ws=w['Contaminant_Reference_List'];names=[]
for row in ws.iter_rows(min_row=5):
 if isinstance(row[0].value,(int,float)) and row[1].value:
  t=str(row[1].value).strip();names.append(dict(reference_row=row[1].row,reference_number=int(row[0].value),reference_taxon=t,in_268=t in x.index,in_41=t in x41.index,in_39=t in x39.index,applied_version_verified=False,contamination_status='Not adjudicated by this inventory'))
refs=pd.DataFrame(names);assert len(refs)==178 and refs.reference_taxon.is_unique and not refs[['in_268','in_41','in_39']].any().any()
refs.to_csv(D/'Table_S2C_1_recovered_reference_names.csv',index=False)
seq=pd.read_csv(I/'Original_seq_clean.csv');junk={'No','No reads','no contaminants','no microbes','no reads','no reads without contaminants','no sample'}
z=seq[~seq.Taxon.isin(junk)&seq.Taxon.notna()&~seq.Taxon.str.lower().str.startswith('no ')].copy();z=z[z.Fraction=='PreSort'].copy();z['SampleID']=z.Patient+'_'+z.Organ
pivot=z.pivot_table(index='Taxon',columns='SampleID',values='Pct',aggfunc='mean').fillna(0)
v=[c for c in pivot if 'Vagina' in c];e=[c for c in pivot if 'Endometrium' in c];c=[c for c in pivot if 'Cervix' in c]
pv=(pivot[v]>0).mean(axis=1);pe=(pivot[e]>0).mean(axis=1)
old=pd.read_csv(I/'Original_otu_clean.csv',index_col=0)
pd.DataFrame(dict(Taxon=pivot.index,vaginal_prevalence=pv,endometrial_prevalence=pe,body_site_prevalence_ratio=pe/(pe+pv+1e-6),in_original_otu_clean=pivot.index.isin(old.index))).to_csv(D/'Historical_body_site_ratio_NOT_contamination_probability.csv',index=False)
script=(I/'Original_sfig16_script.py').read_text();tree=ast.parse(script)
imports=[ast.unparse(n) for n in ast.walk(tree) if isinstance(n,(ast.Import,ast.ImportFrom))]
findings=dict(original_script_sha256=sha(I/'Original_sfig16_script.py'),imports=imports,official_decontam_or_sourcetracker_imported=any('sourcetracker' in s.lower() or 'decontam' in s.lower() for s in imports),
 site_ratio_uses_negative_controls=False,site_ratio_uses_DNA_concentrations=False,historical_pivot_taxa=len(pivot),historical_PRE_body_site_columns=dict(vaginal=len(v),cervical=len(c),endometrial=len(e)),
 model_uses_pooled_fixed_source_profiles='t = S.sum(axis=0)' in script,model_unknown_profile_uniform='un = np.ones(n_taxa)' in script,model_resamples_500_synthetic_assignments='np.random.multinomial(500, sink / total)' in script,model_selects_restart_by_highest_likelihood='if ll > best_ll:' in script,
 interpretation='Custom body-site ratio and source-mixture routine; neither their package equivalence nor their use as independent contamination validation is established. The historical pivot collapses patient/site rows with arithmetic mean, including preparation distinctions when present. No model source proportion is used in the revision. No custom Gibbs reanalysis is claimed.')
assert not findings['official_decontam_or_sourcetracker_imported'] and all(findings[k] for k in ['model_uses_pooled_fixed_source_profiles','model_unknown_profile_uniform','model_resamples_500_synthetic_assignments','model_selects_restart_by_highest_likelihood'])
(D/'Historical_custom_script_audit.json').write_text(json.dumps(findings,indent=2)+'\n')
summary=dict(status='PASS_FOR_NUMERICAL_RECONSTRUCTION',source_libraries=149,exported_taxa=268,reconstructed_selection_taxa=41,submitted_baseline_taxa=39,candidate_reference_names=178,candidate_reference_members_in_268_41_39=[0,0,0],
 reconstructed_rule='At least 0.5 percent of 268-exported-taxon counts in at least two libraries in any one original-code site/fraction stratum; all counts in selected rows retained, with no per-library zeroing.',
 source_and_workbook_denominators_reproduce_41_membership=True,aggregate_exported_counts_discarded_268_to_41_pct=float(100*(1-x41.sum().sum()/x.sum().sum())),
 maximum_library_exported_counts_discarded_268_to_41_pct=float(q.discarded_268_to_41_pct.max()),
 limitations='Reconstruction does not prove historical command execution, applied reference version, true contaminant status or upstream sample identity. The 268 exported assignments are not the complete raw FASTQ community. Taxa/library observations are not independent participants.',
 input_hashes={f.name:sha(f) for f in sorted(I.iterdir()) if f.is_file() and f.name!='Input_origins.json'})
(O/'Source_filter_analysis_provenance.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
