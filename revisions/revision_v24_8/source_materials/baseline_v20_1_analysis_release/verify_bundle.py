#!/usr/bin/env python3
"""Portable, read-only integrity and downstream arithmetic checks.

Python standard library only. Does not recreate raw processing or validate DAPI.
"""
from pathlib import Path
from collections import defaultdict
import csv,json,hashlib,math,statistics

ROOT=Path(__file__).resolve().parent
B=ROOT/'bundle'
def rows(p):
    with p.open(newline='') as f:return list(csv.DictReader(f))
def check(ok,description):
    if not ok:raise ValueError(description)
def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for z in iter(lambda:f.read(1048576),b''):h.update(z)
    return h.hexdigest()
def run():
    hashes=json.loads((ROOT/'SHA256SUMS.json').read_text())
    for name,want in hashes.items():
        p=ROOT/name
        check(p.is_file(),f'Missing file: {name}')
        check(digest(p)==want,f'Changed file: {name}')
    metadata=rows(B/'source_inputs/source_metadata_original_codes.csv')
    check(len(metadata)==149,'Expected 149 retained libraries')
    check(len({x['Participant'] for x in metadata})==15,'Expected 15 sequenced women')
    check(len({x['SpecimenKey'] for x in metadata})==40,'Expected 40 recorded preparations')
    counts=rows(B/'taxon_analysis/data/baseline_39_source_counts.csv')
    first=next(iter(counts[0]));codes=list(counts[0])[1:]
    check(len(counts)==39 and len(codes)==149,'Expected 39 by 149 count matrix')
    check(codes==[x['SourceCode'] for x in metadata],'Library order differs from metadata')
    totals={c:sum(float(r[c]) for r in counts) for c in codes}
    check(all(x>0 for x in totals.values()),'Empty library')
    ra={r[first]:{c:float(r[c])/totals[c] for c in codes} for r in counts}
    pairs=rows(B/'taxon_analysis/results/Matched_IgA_IgM_preparation_manifest.csv')
    saved=rows(B/'taxon_analysis/results/Site_specific_matched_IgA_IgM_representation_scores.csv')
    saved=[x for x in saved if x['taxon_arm']=='baseline_39' and x['rank']=='provisional_species' and x['scope']=='All']
    check(len(saved)==117,'Expected 117 primary site scores')
    for row in saved:
        by_woman=defaultdict(list)
        for p in pairs:
            if p['Location']==row['Location']:by_woman[p['Participant']].append(p)
        A=[];M=[];tax=row['Taxon']
        for pp in by_woman.values():
            A.append(statistics.mean(ra[tax][p['IgA_SourceCode']] for p in pp))
            M.append(statistics.mean(ra[tax][p['IgM_SourceCode']] for p in pp))
        a,m=statistics.mean(A),statistics.mean(M)
        if a+m==0:check(row['Score_ratio_of_matched_woman_means']=='NA','Zero denominator must remain undefined')
        else:check(math.isclose(a/(a+m),float(row['Score_ratio_of_matched_woman_means']),abs_tol=1e-12),'Score differs: '+str((tax,row['Location'])))
        check(sum(a+m>0 for a,m in zip(A,M))==int(row['women_detected_either']),'Detection count differs')
        check(len(A)==int(row['n_matched_women']),'Matched woman count differs')
    gates=rows(B/'facs/data/Table_S2F_6_original_DAPI_region_counts.csv')
    for r in gates:
        parent=int(r['P1_events'])
        for count,pct in [('DAPI_negative_events','DAPI_negative_percent_of_P1'),('DAPI_positive_events','DAPI_positive_percent_of_P1')]:
            check(math.isclose(100*int(r[count])/parent,float(r[pct]),abs_tol=1e-12),'Historical gate arithmetic differs')
    mock=rows(B/'validation/tables/Table_S2V_3_mock_observations_and_summaries.csv')
    for r in mock:
        z=[float(r[f'run_{i}_percent']) for i in range(1,4)]
        check(math.isclose(statistics.mean(z),float(r['mean_percent']),abs_tol=1e-12),'Mock mean differs')
        check(math.isclose(statistics.stdev(z),float(r['sample_SD_percent']),abs_tol=1e-12),'Mock sample SD differs')
    return dict(status='PASS',payload_hashes_checked=len(hashes),retained_libraries=149,sequenced_women=15,recorded_preparations=40,primary_site_scores_recalculated=117,historical_DAPI_region_rows_checked=len(gates),mock_member_summaries_recalculated=len(mock),scope='Integrity and selected downstream arithmetic only; not raw-read reproduction, classifier accuracy, DAPI agreement or bootstrap resampling.')
if __name__=='__main__':
    print(json.dumps(run(),indent=2))
