#!/usr/bin/env python3
"""Run the joint Q15 analysis into a new directory outside the candidate."""
from pathlib import Path
import argparse,subprocess,sys
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--candidate',type=Path,default=Path(__file__).resolve().parents[1])
p.add_argument('--output',type=Path,required=True)
p.add_argument('--rscript',default='Rscript')
p.add_argument('--verify',action='store_true')
a=p.parse_args();root=a.candidate.resolve();out=a.output.resolve()
if root==out or root in out.parents:p.error('Choose output outside the unchanged candidate')
if out.exists():p.error('Choose a new output directory; previous outputs are preserved')
subprocess.run([a.rscript,'--vanilla',str(root/'analysis_scripts/run_Q15_joint_factorial_2026_11_05.R'),str(root/'bundle/joint_factorial_q15'),str(out),'9999','1000'],check=True)
if a.verify:subprocess.run([sys.executable,str(root/'analysis_scripts/verify_Q15_joint_factorial_2026_11_05.py'),str(out),'--report-out',str(out.parent/(out.name+'_verification.json'))],check=True)
