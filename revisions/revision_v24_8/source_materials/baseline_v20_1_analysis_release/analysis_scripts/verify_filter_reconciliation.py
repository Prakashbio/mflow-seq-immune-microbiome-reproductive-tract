#!/usr/bin/env python3
"""Read-only reconstruction of the historical filter reconciliation.

Uses the source-coded matrices and archived tracking/group snapshots, without
claiming to authenticate the original workbook-generation command.
Python standard library only. Writes no files.
"""
from pathlib import Path
from collections import defaultdict
from statistics import mean
import argparse,csv,json,math
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--candidate',type=Path,default=Path(__file__).resolve().parents[1])
a=p.parse_args();root=a.candidate.resolve();b=root/'bundle';checks=0
def rows(p):
 with p.open(newline='') as f:return list(csv.DictReader(f))
def ck(ok,label):
 global checks
 if not ok:raise AssertionError(label)
 checks+=1
mat={};codes=None
for n in [268,41,39]:
 rr=rows(b/'source_filter_inputs'/('Counts_'+str(n)+'_original_codes.csv'))
 first=list(rr[0])[0];cc=list(rr[0])[1:]
 if codes is None:codes=cc
 ck(cc==codes and len(cc)==149 and len(rr)==n,'matrix dimensions/order')
 mat[n]={r[first]:{c:int(r[c]) for c in cc} for r in rr}
mapping=rows(b/'filter_reconciliation/Source_and_legacy_group_mapping.csv')
ck({r['SourceCode'] for r in mapping}==set(codes),'149 source-code mappings')
groups={k:{r['SourceCode']:r[k] for r in mapping} for k in ['source_group','legacy_group']}
den={n:{c:sum(mat[n][t][c] for t in mat[n]) for c in codes} for n in [268,41]}
def passing(n,t,g):
 counts=defaultdict(int)
 for c in codes:
  if mat[n][t][c]*200>=den[n][c]:counts[groups[g][c]]+=1
 return sorted(k for k,v in counts.items() if v>=2)
selected={t for t in mat[268] if passing(268,t,'source_group')}
ck(selected==set(mat[41]),'exact historical 41 membership from 268')
tracking=rows(b/'filter_reconciliation/Workbook_all_41_group_reconciliation.csv')
ck(len(tracking)==41,'41 tracking rows')
for r in tracking:
 t=r['taxon'];legacy=passing(41,t,'legacy_group')
 printed=[] if not r['printed_groups'] else sorted(r['printed_groups'].split('; '))
 ck(legacy==printed and len(legacy)==int(r['printed_groups_count']) and bool(legacy)==(r['printed_per_group_pass']=='Yes'),'legacy flags '+t)
 ck(passing(268,t,'source_group')==r['original_source_code_groups_268'].split('; '),'source groups '+t)
six=rows(b/'filter_reconciliation/Workbook_six_filter_flags_reconciled.csv')
ck(len(six)==6,'six No/retained cases')
for r in six:ck(r['printed_per_group_pass']=='No' and r['retained_in_39']=='True' and r['taxon'] in mat[39] and bool(passing(268,r['taxon'],'source_group')),'retained passing '+r['taxon'])
numeric=rows(b/'filter_reconciliation/Workbook_all_287_global_numeric_checks.csv')
ck(len(numeric)==287,'287 archived global cells')
for r in numeric:
 t=r['taxon'];v=[mat[41][t][c] for c in codes];ra=[100*mat[41][t][c]/den[41][c] for c in codes]
 calc={'total_counts':sum(v),'detected_libraries':sum(x>0 for x in v),'prevalence_percent':round(100*sum(x>0 for x in v)/149,1),'mean_percent':round(mean(ra),4),'max_percent':round(max(ra),2),'n_ge_0_5_percent':sum(mat[41][t][c]*200>=den[41][c] for c in codes),'n_ge_1_percent':sum(mat[41][t][c]*100>=den[41][c] for c in codes)}[r['field']]
 ck(math.isclose(calc,float(r['recorded']),abs_tol=1e-8,rel_tol=0),'global cell '+r['cell'])
print(json.dumps({'status':'PASS','checks_passed':checks,'starting_assignments':268,'reconstructed_intermediate_assignments':41,'unchanged_baseline_assignments':39,'legacy_group_rows_reproduced':41,'global_numeric_cells_reproduced':287,'six_retained_taxa_pass_original_code_rule':True,'files_written':0,'scope':'Numerical reconstruction from supplied matrices and tracking snapshots; not proof of original execution or contamination status.'},indent=2))
