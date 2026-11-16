"""Separate label correspondence, matching and weighting effects on scores.
No original data or labels are edited. Count-vector correspondence is not a
claim that the women named by different headers are the same person.
"""
from pathlib import Path
import argparse,hashlib,json
import numpy as np
import pandas as pd
R=Path(__file__).resolve().parents[1]
a=argparse.ArgumentParser();a.add_argument('--package-root',type=Path);a.add_argument('--output',type=Path);args=a.parse_args()
O=args.output or R/'taxon_integrated_working_v3_2026-10-28';D=O/'taxon_audit';D.mkdir(parents=True,exist_ok=True)
src=R/'03_analysis/Taxon_scores_PRE_review_2026-10-28'
if args.package_root:
 p=args.package_root/'taxon_inputs';oldpath=p/'historical_otu_table.csv';mdpath=p/'historical_sample_metadata.csv';x=pd.read_csv(p/'baseline_39_source_counts.csv',index_col=0);m=pd.read_csv(p/'source_metadata_with_vaginal_category.csv');pairs=pd.read_csv(args.package_root/'taxon_statistics/Matched_IgA_IgM_preparation_manifest.csv');cross=pd.read_csv(p/'Independent_149_original_report_vector_crosswalk.csv')
else:
 p=R.parent/'mFLOW-seq_Menuscript/Bioinfomatics analyse by KASUN /Third claude/mFLOW-seq_Menuscript/Scripts/final_scripts_1803_new/data'
 oldpath=p/'otu_table.csv';mdpath=p/'sample_metadata.csv';x=pd.read_csv(src/'data/baseline_39_source_counts.csv',index_col=0);m=pd.read_csv(src/'data/source_metadata_with_vaginal_category.csv');pairs=pd.read_csv(src/'results/Matched_IgA_IgM_preparation_manifest.csv');cross=pd.read_csv(R/'02_evidence/Project_folder_audit_2026-10-27/Independent_149_original_report_vector_crosswalk.csv')
old=pd.read_csv(oldpath,index_col=0);om=pd.read_csv(mdpath).set_index('SampleID').loc[old.columns]
assert old.shape==x.shape==(39,149) and np.array_equal(old.to_numpy(),x.to_numpy())
assert list(m.SourceCode)==list(x.columns)
for row in cross.itertuples():assert np.array_equal(x[row.source_code].to_numpy(),old[row.original_table_heading].to_numpy())
pr=old.div(old.sum(),axis=1);sr=x.div(x.sum(),axis=1)
ratio=lambda aa,bb:np.divide(aa,aa+bb,out=np.full(len(aa),np.nan),where=(aa+bb)>0)
global_old=ratio(pr.loc[:,om.Antibody.eq('A')].mean(axis=1).to_numpy(),pr.loc[:,om.Antibody.eq('M')].mean(axis=1).to_numpy())
global_source=ratio(sr[m.loc[m.Fraction.eq('IgA'),'SourceCode']].mean(axis=1).to_numpy(),sr[m.loc[m.Fraction.eq('IgM'),'SourceCode']].mean(axis=1).to_numpy())
out=[];contributions=[]
for site in ['VAGINA','CERVIX','ENDOMETRIUM']:
 oa=pr.loc[:,om.Location.eq(site)&om.Antibody.eq('A')].mean(axis=1).to_numpy();ob=pr.loc[:,om.Location.eq(site)&om.Antibody.eq('M')].mean(axis=1).to_numpy()
 sa=sr[m.loc[m.Location.eq(site)&m.Fraction.eq('IgA'),'SourceCode']].mean(axis=1).to_numpy();sb=sr[m.loc[m.Location.eq(site)&m.Fraction.eq('IgM'),'SourceCode']].mean(axis=1).to_numpy()
 z=pairs[pairs.Location==site];pa=sr[z.IgA_SourceCode].mean(axis=1).to_numpy();pb=sr[z.IgM_SourceCode].mean(axis=1).to_numpy()
 A=[];M=[]
 for woman,zz in z.groupby('Participant'):
  A.append(sr[zz.IgA_SourceCode].mean(axis=1));M.append(sr[zz.IgM_SourceCode].mean(axis=1))
 ma=np.mean(A,axis=0);mb=np.mean(M,axis=0)
 oldsite,srcsite,paired,weighted=[ratio(a,b) for a,b in [(oa,ob),(sa,sb),(pa,pb),(ma,mb)]]
 for i,tax in enumerate(x.index):
  out.append(dict(Taxon=tax,Location=site,historical_global_pooled_score=global_old[i],named_source_global_pooled_score=global_source[i],historical_site_unmatched_score=oldsite[i],named_source_site_unmatched_score=srcsite[i],named_source_site_matched_preparation_weighted_score=paired[i],named_source_site_matched_woman_weighted_score=weighted[i],delta_label_correspondence=srcsite[i]-oldsite[i],delta_pair_eligibility=paired[i]-srcsite[i],delta_woman_weighting=weighted[i]-paired[i],matched_preparations=len(z),matched_women=z.Participant.nunique()))
bridge=pd.DataFrame(out);bridge.to_csv(D/'Score_definition_and_correspondence_bridge.csv',index=False)
saved=pd.read_csv((args.package_root/'taxon_statistics' if args.package_root else src/'results')/'Site_specific_matched_IgA_IgM_representation_scores.csv')
saved=saved[(saved.taxon_arm=='baseline_39')&(saved['rank']=='provisional_species')&(saved.scope=='All')]
both=bridge.merge(saved,on=['Taxon','Location'],validate='one_to_one');assert len(both)==117
assert np.allclose(both.named_source_site_matched_woman_weighted_score,both.Score_ratio_of_matched_woman_means,atol=1e-12,equal_nan=True)
# Trace the original count vectors carrying the largest Sneathia abundances.
for code in sr.columns:
 c=cross[cross.source_code==code].iloc[0];meta=m[m.SourceCode==code].iloc[0];hist=om.loc[c.original_table_heading]
 contributions.append(dict(source_code=code,source_participant=int(meta.Participant),source_site=meta.Location,source_fraction=meta.Fraction,historical_table_heading=c.original_table_heading,historical_metadata_participant=str(hist.Participant) if 'Participant' in hist else '',historical_site=hist.Location,historical_fraction=hist.Antibody,Sneathia_relative_abundance=sr.loc['Sneathia vaginalis',code],Sneathia_retained_counts=int(x.loc['Sneathia vaginalis',code]),retained_library_total=int(x[code].sum()),report=c.report,report_sha256=c.report_sha256,interpretation='Count-vector lineage only; original labels remain unchanged; differing participant labels are not equated.'))
pd.DataFrame(contributions).sort_values('Sneathia_relative_abundance',ascending=False).to_csv(D/'Sneathia_library_count_vector_trace.csv',index=False)
summary=dict(status='PASS',libraries=149,source_vs_historical_value_matrix_identical=True,all_149_crosswalk_vectors_checked=True,site_taxon_rows=117,final_scores_match_independent_reconstruction=True,original_files_unchanged=True,limitation='Historical headings and named source codes allocate identical count vectors to different factor labels. This is a provenance hold, not merely a model adjustment or permission to rename participants. Site stratification, pairing and woman weighting are additional, separately quantified changes. Upstream identity and final author adjudication remain unresolved.')
(D/'Verification.json').write_text(json.dumps(summary,indent=2)+'\n')
(D/'READ_FIRST.txt').write_text('''TAXON SCORE PROVENANCE BRIDGE —28 October2026

The same39×149 numeric matrix is present under historical derived headings
and named original source codes. The score bridge holds the counts fixed and
separates site/fraction-label allocation, IgA/IgM pairing eligibility and equal
woman weighting. Every final site score matches the checked taxon analysis.

These comparisons must NOT be described as a change caused solely by paired
statistical analysis. The original-source/derived-heading correspondence still
needs author adjudication. A vector association does not mean participant104
and90 are the same woman. No original names or input files are changed.

For example, the historical pooled Sneathia vaginalis score is0.999482.
The named-source, matched-woman site scores are0.406411(vagina),0.268873(cervix)
and0.010993(endometrium). The bridge provides intermediate calculations, and
the library trace identifies every source vector and its original report hash.
These estimates remain provisional and do not validate taxonomic accuracy or
establish that a biological conclusion is correct/incorrect before provenance
is resolved. They are not per-cell coating measurements.
''')
print(json.dumps(summary,indent=2))
