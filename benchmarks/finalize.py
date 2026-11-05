"""Final controlled optimization comparison, artifact checks and report creation."""
import json
import subprocess
import sys
import time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]


def run(parts,log=None):
    t=time.perf_counter()
    if log:
        with (ROOT/log).open('w',encoding='utf-8') as f:
            subprocess.run([sys.executable,*parts],cwd=ROOT,stdout=f,stderr=subprocess.STDOUT,check=True)
    else: subprocess.run([sys.executable,*parts],cwd=ROOT,check=True)
    return time.perf_counter()-t


print('FINAL correctness and output-copy optimization',flush=True)
run(['tests/test_correctness.py'],'tests/test_results.txt')
run(['benchmarks/benchmark.py','--project','1','--scope','screening','--methods','radix-parallel','--repeats','3'])
run(['benchmarks/benchmark.py','--project','1','--scope','full','--methods','radix-parallel','radix-jit-threaded','quicksort','--repeats','3','--keep-method','radix-parallel'])
print('FINAL strict reference sampling',flush=True)
run(['src/rescale.py','--policy','strict','--output',str(ROOT/'results/rescale/query.strict.bed')],'results/rescale/strict_run.log')
print('FINAL figures and patient analysis',flush=True)
times=[]
for project in [2,3,4]:
    elapsed=run(['src/make_figures.py','--project',str(project)],f'benchmarks/figures_p{project}.log')
    times.append(dict(project=project,figure_process_seconds=elapsed,includes_statistics=project==4,repeats=1,input='previously generated output table/matrix',separate_from_algorithm_ranking=True))
(ROOT/'benchmarks/figure_times.json').write_text(json.dumps(times,indent=2))
sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from vplot import write_matrix
write_matrix(ROOT/'results/vplot/vplot.long.tsv',np.load(ROOT/'results/vplot/vplot.matrix.npy'),-500,True)
print('FINAL independent verification',flush=True)
run(['tests/validate_outputs.py'],'tests/output_validation.txt')
run(['benchmarks/report.py'])
print('FINAL reports generated',flush=True)
