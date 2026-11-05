"""Recorded sequential campaign for the remaining alternatives."""
import subprocess
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
jobs=[
    ['prepare_data.py'],
    ['tests/test_correctness.py'],
    ['benchmarks/benchmark.py','--project','1','--scope','screening','--methods','radix-threaded','radix-jit','radix-jit-threaded','--repeats','3'],
    ['benchmarks/benchmark.py','--project','1','--scope','full','--methods','radix-threaded','radix-jit','radix-jit-threaded','--repeats','3','--keep-method','radix-jit-threaded'],
]
for project,keep in [(2,'compiled'),(3,'compiled'),(4,'bytes')]:
    jobs.append(['benchmarks/benchmark.py','--project',str(project),'--scope','screening','--repeats','3'])
    jobs.append(['benchmarks/benchmark.py','--project',str(project),'--scope','full','--repeats','3','--keep-method',keep])
jobs.append(['benchmarks/jit_cost.py'])
for job in jobs:
    print('CAMPAIGN', ' '.join(job),flush=True)
    subprocess.run([sys.executable,*job],cwd=ROOT,check=True)
