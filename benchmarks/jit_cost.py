"""Measure tiny first calls with an explicitly empty, dedicated Numba cache."""
import json
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]


def worker():
    sys.path.insert(0,str(ROOT/'src'))
    import numpy as np
    from kernels import bed_records,bed_offsets_lengths,gather,stable_ties,exact_select,vplot_rows
    from sort_kernels import radix_indices
    from parallel_writer import copy_rows
    b=np.frombuffer(bytearray(b'chr1\t1\t3\nchr1\t4\t8\n'),np.uint8)
    n=np.array([[99,104,114,49]],np.uint8); lens=np.array([4],np.int64)
    calls=[('BED parser',lambda:bed_records(b,n,lens)),('length parser',lambda:bed_offsets_lengths(b)),('tie ordering',lambda:stable_ties(np.array([1,0],np.int64),np.array([5,5],np.uint64))),('exact sampling',lambda:exact_select(np.array([1,1],np.int64),np.array([0,2],np.int64),np.array([0,1],np.int64),2026)),('V-plot',lambda:vplot_rows(np.frombuffer(b'2\tchr1\t0\t10\tm\t1\t+\tchr1\t0\t8\tr\t8\t+\n',np.uint8),-500,499,1000,True))]
    rec=np.array([[0,1,3,0,9],[0,4,8,9,9]],np.uint64)
    calls.append(('BED output copy',lambda:gather(b,rec[:,3],rec[:,4],np.array([1,0],np.int64),0,np.empty(1024,np.uint8))))
    calls.append(('Parallel BED copy',lambda:copy_rows(b,rec[:,3],rec[:,4],np.array([1,0],np.int64),np.array([0,9,18],np.int64),np.empty(18,np.uint8))))
    calls.append(('Explicit radix sort',lambda:radix_indices(np.array([99,2,99,0],np.uint64))))
    records=[]
    for name,call in calls:
        t=time.perf_counter(); call(); first=time.perf_counter()-t
        t=time.perf_counter(); call(); second=time.perf_counter()-t
        records.append(dict(kernel=name,first_call_seconds=first,in_process_second_call_seconds=second))
    print(json.dumps(records))


if __name__=='__main__':
    if '--worker' in sys.argv: worker()
    else:
        with tempfile.TemporaryDirectory(prefix='jit-cost-',dir=ROOT/'tmp') as temp:
            env=dict(os.environ,NUMBA_CACHE_DIR=temp,NUMBA_NUM_THREADS='1')
            runs=[]
            for mode in ['empty_disk_cache','populated_disk_cache']:
                t=time.perf_counter()
                result=subprocess.run([sys.executable,__file__,'--worker'],env=env,capture_output=True,text=True,check=True)
                runs.append(dict(mode=mode,process_seconds=time.perf_counter()-t,kernels=json.loads(result.stdout)))
            (ROOT/'benchmarks/jit_cost.json').write_text(json.dumps(runs,indent=2))
            print(json.dumps(runs,indent=2))
