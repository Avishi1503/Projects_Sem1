"""One full-data empty-JIT-cache check of each selected Numba-based winner."""
import json
import os
import tempfile
from pathlib import Path
from benchmark import ROOT,run_job

winners=json.loads((ROOT/'benchmarks/winners.json').read_text())
rows=[]
for project in [1,2,3]:
    with tempfile.TemporaryDirectory(prefix=f'cold-p{project}-',dir=ROOT/'tmp') as temporary:
        cache=Path(temporary).resolve()
        assert cache.is_relative_to((ROOT/'tmp').resolve())
        previous=os.environ.get('NUMBA_CACHE_DIR')
        os.environ['NUMBA_CACHE_DIR']=str(cache)
        try:
            print('EMPTY JIT CACHE',project,winners[str(project)]['method'],flush=True)
            row=run_job(project,winners[str(project)]['method'],'full',101,False)
        finally:
            if previous is None: os.environ.pop('NUMBA_CACHE_DIR',None)
            else: os.environ['NUMBA_CACHE_DIR']=previous
        if row['exit_code']: raise RuntimeError(row['error'])
        row['numba_cache']='empty at process launch'
        row['repeats']=1
        rows.append(row)
        print(json.dumps({k:row[k] for k in ['project','method','wall_seconds','output_sha256']}),flush=True)
(ROOT/'benchmarks/cold_start_full.json').write_text(json.dumps(rows,indent=2))
