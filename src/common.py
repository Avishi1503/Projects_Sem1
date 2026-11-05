import contextlib
import gzip
import json
import sys
import time
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def json_write(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(obj,indent=2,allow_nan=False),encoding='utf-8')


@contextlib.contextmanager
def open_input(path):
    if str(path) == '-':
        yield sys.stdin.buffer
    else:
        with (gzip.open(path,'rb') if str(path).endswith('.gz') else open(path,'rb')) as f:
            yield f


def read_bytes(paths, limit=0):
    """Read gzip (including concatenated members) or plain BED, preserving lines."""
    out = bytearray()
    remaining = limit
    for path in paths:
        with open_input(path) as f:
            if limit:
                for line in f:
                    out.extend(line)
                    remaining -= 1
                    if remaining == 0:
                        break
            else:
                for chunk in iter(lambda:f.read(8*1024*1024),b''):
                    out.extend(chunk)
        if out and out[-1] != 10:
            out.extend(b'\n')
        if limit and remaining == 0:
            break
    return np.frombuffer(out,dtype=np.uint8)


def selection(path):
    names = [x.strip().encode() for x in Path(path).read_text().splitlines() if x.strip() and not x.startswith('#')]
    if not names or len(set(names)) != len(names):
        raise ValueError('Selection must contain unique chromosome names')
    if any(b'\t' in s or b' ' in s for s in names):
        raise ValueError('Selection must contain one chromosome per line')
    array = np.zeros((len(names),max(map(len,names))),np.uint8)
    for i,n in enumerate(names):
        array[i,:len(n)] = np.frombuffer(n,np.uint8)
    return names,array,np.array(list(map(len,names)),np.int64)


def write_order(path,b,offsets,lengths,order):
    from kernels import gather
    path = Path(path)
    path.parent.mkdir(parents=True,exist_ok=True)
    buf = np.empty(8*1024*1024,np.uint8)
    begin = 0
    with path.open('wb') as f:
        while begin < len(order):
            begin,used = gather(b,offsets,lengths,order,begin,buf)
            f.write(memoryview(buf[:used]))


def emit_metrics(args, info, started):
    info['pipeline_seconds'] = time.perf_counter()-started
    if getattr(args,'metrics',None):
        json_write(args.metrics,info)
    print(json.dumps(info),flush=True)
