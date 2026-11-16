"""Recompute v25 supplied-summary arithmetic using Python's standard library.

Usage: python3 recompute_v25_sequencing_summaries.py --input data --output output
This does NOT rerun Kraken2, Bracken, simulation generation or clinical models.
"""
from pathlib import Path
import argparse, csv, json, statistics
from collections import defaultdict

def main():
 ap=argparse.ArgumentParser(description=__doc__)
 ap.add_argument('--input',type=Path,default=Path(__file__).resolve().parent.parent/'data')
 ap.add_argument('--output',type=Path,default=Path('v25_summary_check_output'))
 a=ap.parse_args();a.output.mkdir(parents=True,exist_ok=True)
 read=lambda n:list(csv.DictReader((a.input/n).open(newline='')))
 checks=[]
 def ck(label,value):checks.append({'check':label,'passed':bool(value)})
 def val(r,k):return float(r[k])
 def stats(v):return dict(min=min(v),median=statistics.median(v),max=max(v))
 q=read('Supplementary_QC_source_values.csv')
 ck('149 distinct supplied library codes',len(q)==149==len({r['Sample ID'] for r in q}))
 for r in q:
  v=lambda k:val(r,k)
  ck(r['Sample ID']+' received balance',v('Reads received')==v('Reads classified (Kraken2)')+v('Reads unclassified'))
  ck(r['Sample ID']+' classified ranks balance',v('Reads classified (Kraken2)')==sum(v(k) for k in ['Direct species-level reads','Direct genus-level reads','Reads above genus level']))
  ck(r['Sample ID']+' Bracken stage balance',v('Bracken: total species-level reads')==v('Bracken: directly assigned reads')+v('Bracken: redistributed reads'))
  for p,n,d in [('Classified (%)','Reads classified (Kraken2)','Reads received'),('Unclassified (%)','Reads unclassified','Reads received'),('Redistributed share of estimate (%)','Bracken: redistributed reads','Bracken: total species-level reads')]:
   ck(r['Sample ID']+' '+p,abs(v(p)-100*v(n)/v(d))<=.005001)
 single=read('benchmark_single_species.csv');ck('Nine single-species sources',len(single)==9)
 single_out=[]
 for r in single:
  n=val(r,'total_reads');correct=val(r,'kraken2_correct_species_reads');wrong=val(r,'kraken2_wrong_species_reads')
  cats=['kraken2_correct_species_reads','kraken2_wrong_species_reads','kraken2_genus_only_reads','kraken2_above_genus_reads','unclassified_reads']
  ck(r['species']+' Kraken partition',sum(val(r,k) for k in cats)==n==100000)
  bc=val(r,'bracken_correct_reads');bw=val(r,'bracken_wrong_reads')
  ck(r['species']+' Bracken correct denominator',abs(100*bc/(bc+bw)-val(r,'bracken_correct_pct'))<=.005001)
  single_out.append(dict(species=r['species'],correct_Kraken_percent_input=100*correct/n,
   incorrect_Kraken_percent_input=100*wrong/n,incorrect_Kraken_percent_tabulated_species=100*wrong/(correct+wrong),
   correct_Bracken_percent_input=100*bc/n,correct_Bracken_percent_listed_estimates=100*bc/(bc+bw)))
 equal=read('benchmark_mock_community.csv');ck('Nine equal-mixture source taxa',len(equal)==9)
 m=read('mixture_results.csv');groups=defaultdict(list)
 for r in m:groups[r['sample']].append(r)
 ck('56 rows,28 pair mixtures',len(m)==56 and len(groups)==28)
 mix_out=[];definition_flags=[]
 for sample,rr in sorted(groups.items()):
  ck(sample+' two taxa sum to100 input',len(rr)==2 and sum(val(r,'expected_pct') for r in rr)==100)
  unequal=len({r['expected_pct'] for r in rr})==2
  preserved=([r['species'] for r in sorted(rr,key=lambda r:val(r,'expected_pct'))]==[r['species'] for r in sorted(rr,key=lambda r:val(r,'bracken_pct'))]) if unequal else None
  mix_out.append(dict(sample=sample,pair=rr[0]['pair'],unequal_input=unequal,rank_order_preserved=preserved,
   complement_not_assigned_to_two_input_taxa=100-sum(val(r,'bracken_pct') for r in rr)))
  for r in rr:
   ck(sample+' signed source abs_error_pp '+r['species'],abs(val(r,'abs_error_pp')-(val(r,'bracken_pct')-val(r,'expected_pct')))<=.005001)
   delta=val(r,'bracken_est_reads')-val(r,'bracken_redistributed_reads')-val(r,'kraken2_direct_reads')
   if delta:definition_flags.append(dict(sample=sample,species=r['species'],implied_clade_minus_tabulated_direct=delta))
 for name,rows in [('Single_species_recomputed.csv',single_out),('Pair_mixture_recomputed.csv',mix_out),('Count_definition_flags.csv',definition_flags),('Arithmetic_checks.csv',checks)]:
  with (a.output/name).open('w',newline='') as f:
   w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 unequal=[r for r in mix_out if r['unequal_input']]
 report=dict(status='PASS' if all(r['passed'] for r in checks) else 'ARITHMETIC_FINDING',checks=len(checks),
  QC={k:stats([val(r,k) for r in q]) for k in ['Reads received','Reads classified (Kraken2)','Classified (%)','Redistributed share of estimate (%)']},
  low_depth=[dict(sample=r['Sample ID'],classified=int(val(r,'Reads classified (Kraken2)'))) for r in q if val(r,'Reads classified (Kraken2)')<10000],
  equal_mixture_listed_percent_sum=sum(val(r,'bracken_est_pct') for r in equal),
  equal_mixture_complement=100-sum(val(r,'bracken_est_pct') for r in equal),
  unequal_mixtures=len(unequal),ranks_preserved=sum(r['rank_order_preserved'] for r in unequal),
  one_percent_observed_range=stats([val(r,'bracken_pct') for r in m if val(r,'expected_pct')==1]),
  distinct_count_definition_flags=len(definition_flags),
  limitations=['These computations audit supplied summaries, not independent assay accuracy.','Count-definition flags require original clade/direct report definitions; source values are not corrected.','No sequencing classifier, simulation generator or clinical model is executed.'])
 (a.output/'Summary.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps(report,indent=2))
 assert report['status']=='PASS'
if __name__=='__main__':main()
