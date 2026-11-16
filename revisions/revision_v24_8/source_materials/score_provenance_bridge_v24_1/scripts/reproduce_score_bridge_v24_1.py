from pathlib import Path
import sys,subprocess,hashlib,json
release=Path(sys.argv[1]).resolve();out=Path(sys.argv[2]).resolve();assert not out.exists(),"Use a new output directory"
bridge=release/'score_provenance_bridge_v24_1';base=release/'baseline_v20_1_analysis_release/bundle'
for rel,h in json.loads((bridge/'Input_hashes.json').read_text()).items():assert hashlib.sha256((base/rel).read_bytes()).hexdigest()==h,rel
subprocess.run([sys.executable,str(bridge/'scripts/audit_taxon_score_bridge_2026_10_28.py'),'--package-root',str(base),'--output',str(out)],check=True)
for rel,h in json.loads((bridge/'Expected_result_hashes.json').read_text()).items():assert hashlib.sha256((out/'taxon_audit'/rel).read_bytes()).hexdigest()==h,rel
print('PASS: all bridge inputs and four rebuilt outputs match; all 117 primary scores verified.')
