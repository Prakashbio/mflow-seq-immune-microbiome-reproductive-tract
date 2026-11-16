"""New review copies from original exports and workbook-confirmed keys.
No original labels/data edited; no Kasun analysis input used.
"""
from pathlib import Path
import os
import shutil,json,hashlib
import numpy as np
import pandas as pd
import openpyxl
R=Path(__file__).resolve().parents[1];O=R/'03_analysis/Taxon_scores_PRE_review_2026-10-28'
for d in ['data','results','documents']:(O/d).mkdir(parents=True,exist_ok=True)
S=R/'03_analysis/Prakash_original_workflow_review_2026-10-27/data'
P=R.parent/'mFLOW-seq_Menuscript/Bioinfomatics analyse by KASUN /Third claude/mFLOW-seq_Menuscript/Scripts/final_scripts_1803_new'
manifest=R/'02_evidence/Source_records_2026-10-20/Source_record_manifest_REVIEW_ONLY.csv'
book=Path(os.environ.get('MFLOW_NATIVE_WORKBOOK', 'Prakash FACS study final version of microbiome sequencing. Improved pipline 07.04.25.xlsx'))
paths=[manifest,book,P/'data/otu_table.csv',P/'data/sample_metadata.csv',P/'data/taxonomy.csv',P/'Results/SupplementaryTableS1_IgACoatingScores.csv']+[S/n for n in ['baseline_39_source_counts.csv','sensitivity_41_source_counts.csv','source_metadata_original_codes.csv','Workbook_preparation_key_evidence.csv']]
hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
for p in paths[6:]:shutil.copy2(p,O/'data'/p.name)
shutil.copy2(P/'data/taxonomy.csv',O/'data/original_taxonomy.csv')
m=pd.read_csv(S/'source_metadata_original_codes.csv');fm=pd.read_csv(manifest);w=openpyxl.load_workbook(book,data_only=True,read_only=False)
cats=[]
for person,g in fm.groupby('Participant'):
 descriptions=set(g.workbook_community_description);assert len(descriptions)==1,(person,descriptions)
 desc=next(iter(descriptions));cells=sorted(set('E'+x.split('!D')[-1] for x in g.workbook_patient_cell))
 for c in cells:assert str(w['General table'][c].value)==desc,(person,c)
 category='VAGINAL_DYSBIOSIS_DESCRIPTION' if 'dysbiosis' in desc.lower() else 'LACTOBACILLUS_DESCRIPTION'
 assert category=='VAGINAL_DYSBIOSIS_DESCRIPTION' or any(x in desc.lower() for x in ['iner','crispatus','jenn','jens'])
 cats.append(dict(Participant=int(person),VaginalCategory=category,Original_workbook_description=desc,Original_description_cells='General table!'+','.join(cells),classification_scope='Original vaginal-community description, not clinical disease or independently validated CST.'))
w.close();cats=pd.DataFrame(cats);cats.to_csv(O/'data/Participant_vaginal_category_from_native_workbook.csv',index=False)
m=m.merge(cats[['Participant','VaginalCategory']],on='Participant',how='left',validate='many_to_one');assert m.VaginalCategory.notna().all() and m.SpecimenKey.nunique()==40
assert m.set_index('SourceCode').loc[['78EA2','78EG2','78EM2','78EP'],'SpecimenKey'].nunique()==1
m.to_csv(O/'data/source_metadata_with_vaginal_category.csv',index=False)
# Reproduce the original score calculations arithmetically without importing
# the historical scripts or adopting their biological interpretation.
old=pd.read_csv(P/'data/otu_table.csv',index_col=0);md=pd.read_csv(P/'data/sample_metadata.csv').set_index('SampleID');rel=old.div(old.sum(),axis=1)
common=rel.columns.intersection(md.index);rel=rel[common];md=md.loc[common];assert len(common)==149
a=rel.loc[:,md.Antibody.eq('A')].mean(axis=1);b=rel.loc[:,md.Antibody.eq('M')].mean(axis=1)
score=a/(a+b).replace(0,np.nan);by=[]
for site in ['VAGINA','CERVIX','ENDOMETRIUM']:
 aa=rel.loc[:,md.Antibody.eq('A')&md.Location.eq(site)].mean(axis=1);bb=rel.loc[:,md.Antibody.eq('M')&md.Location.eq(site)].mean(axis=1);by.append(aa/(aa+bb).replace(0,np.nan))
equal=pd.concat(by,axis=1).mean(axis=1);hist=pd.read_csv(P/'Results/SupplementaryTableS1_IgACoatingScores.csv').set_index('Taxon')
result=pd.DataFrame({'Taxon':score.index,'Original_pooled_ratio_of_group_means':score.values,'Original_equal_available_site_mean':equal.values,'Exported_original_table_score':hist.loc[score.index,'IgA_coating_score'].values})
result['absolute_difference_export_vs_global_formula']=(result.Original_pooled_ratio_of_group_means-result.Exported_original_table_score).abs();result.to_csv(O/'results/Historical_score_formula_reproduction.csv',index=False)
assert result.absolute_difference_export_vs_global_formula.max()<0.00050001
summary=dict(original_source_hashes=hashes,libraries=len(m),women=m.Participant.nunique(),workbook_preparation_keys=m.SpecimenKey.nunique(),category_women=cats.VaginalCategory.value_counts().to_dict(),original_score_reproduction_max_difference=float(result.absolute_difference_export_vs_global_formula.max()),score_provenance='Historical table uses pooled ratio of means; original sfig10 orders taxa using mean available-site scores. These are different estimands. Reproduction does not validate original derived headings or per-cell interpretation.',scope='Separate exploratory revision analysis from named source exports. Original source headings and count values unchanged; laboratory/run linkage and final source-derived label correspondence remain unresolved.')
assert all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in hashes.items())
(O/'Input_provenance.json').write_text(json.dumps(summary,indent=2));print(json.dumps({k:v for k,v in summary.items() if k!='original_source_hashes'},indent=2))
