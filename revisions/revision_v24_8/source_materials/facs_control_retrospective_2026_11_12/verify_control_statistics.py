"""Check reported control arithmetic; counts are visually transcribed, not inferred gates."""
from pathlib import Path
from decimal import Decimal,ROUND_HALF_UP
import csv,hashlib,json
r=Path(__file__).resolve().parent
with (r/'New_HUT093_PDF_population_transcriptions.csv').open(newline='') as f:rows=list(csv.DictReader(f))
for x in rows:
 d=Decimal(x['DAPI_labelled_negative']);m=Decimal(x['IgM_region']);value=100*m/d
 assert value.quantize(Decimal('0.1'),rounding=ROUND_HALF_UP)==Decimal(x['printed_IgM_parent_percent'])
 assert abs(float(value)-float(x['IgM_exact_parent_percent']))<1e-12
 assert int(x['P1'])-int(d)-int(x['DAPI_labelled_positive'])==int(x['P1_outside_both_DAPI_regions'])
 assert int(d)-sum(int(x[k]) for k in ['IgM_region','IgA_region','Double_positive_region','negative_rectangle'])==int(x['DAPI_negative_outside_four_antibody_rectangles'])
 assert hashlib.sha256((r/'source_reports'/x['source_file']).read_bytes()).hexdigest()==x['source_sha256']
with (r/'Historical_vs_retrospective_control_statistics.csv').open(newline='') as f:history=list(csv.DictReader(f))
for x in history:
 for stage in ['historical','retrospective']:
  v=100*int(x[stage+'_IgM_region'])/int(x[stage+'_DAPI_labelled_parent']);assert abs(v-float(x[stage+'_exact_IgM_parent_percent']))<1e-12
print(json.dumps(dict(status='PASS',retrospective_reports=len(rows),historical_comparisons=len(history),source_PDF_hashes_verified=True,clinical_data_changes=0)))
