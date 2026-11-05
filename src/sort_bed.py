"""Project 1: preserve full BED lines, filter and sort in selection-file order."""
import argparse
import gc
import heapq
import tempfile
import time
from pathlib import Path
import numpy as np
from common import ROOT, selection, read_bytes, write_order, emit_metrics, open_input
from kernels import bed_records, stable_ties

METHODS = ['radix','radix-threaded','radix-jit','radix-jit-threaded','radix-parallel','quicksort','heapsort','lexsort','timsort','bucket','external']


def packed_keys(records):
    max_start = int(records[:,1].max(initial=0))
    max_end = int(records[:,2].max(initial=0))
    sb = max(1,max_start.bit_length())
    eb = max(1,max_end.bit_length())
    rb = int(records[:,0].max(initial=0)).bit_length()
    if rb+sb+eb > 64:
        return None
    return (records[:,0] << np.uint64(sb+eb)) | (records[:,1] << np.uint64(eb)) | records[:,2]


def sorted_order(records,method):
    n = len(records)
    if method == 'lexsort':
        return np.lexsort((records[:,2],records[:,1],records[:,0]))
    if method == 'timsort':
        return np.array(sorted(range(n),key=lambda i:(int(records[i,0]),int(records[i,1]),int(records[i,2]))),np.int64)
    if method == 'bucket':
        buckets = [[] for _ in range(int(records[:,0].max(initial=0))+1)]
        for i in range(n):
            buckets[int(records[i,0])].append(i)
        for bucket in buckets:
            bucket.sort(key=lambda i:(int(records[i,1]),int(records[i,2])))
        return np.fromiter((i for bucket in buckets for i in bucket),dtype=np.int64,count=n)
    keys = packed_keys(records)
    if keys is None:
        # No truncation or overflow for wider coordinates/large selections.
        return np.lexsort((records[:,2],records[:,1],records[:,0]))
    if method in ('radix-jit','radix-jit-threaded','radix-parallel'):
        from sort_kernels import radix_indices
        return radix_indices(keys)
    order = np.argsort(keys,kind='stable' if method in ('radix','radix-threaded') else method)
    if method not in ('radix','radix-threaded'):
        stable_ties(order,keys)
    return order


def external_sort(paths,selection_path,output,chunk_rows=250000,limit=0):
    names,_,_ = selection(selection_path)
    ranks = {n:i for i,n in enumerate(names)}
    seen = kept = 0
    def key(line):
        fields = line.split(b'\t',3)
        return ranks[fields[0]],int(fields[1]),int(fields[2])
    output = Path(output)
    output.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='bed-sort-',dir=output.parent) as temporary:
        chunks = []
        lines = []
        def spill():
            lines.sort(key=key)
            path = Path(temporary)/f'{len(chunks):05d}.bed'
            with path.open('wb') as f:
                f.writelines(lines)
            chunks.append(path)
            lines.clear()
        input_rows = 0
        for path in paths:
            with open_input(path) as f:
                for line in f:
                    input_rows += 1
                    if limit and input_rows > limit:
                        break
                    line = line.rstrip(b'\r\n')
                    if not line or line.startswith((b'#',b'track ',b'browser ')):
                        continue
                    fields = line.split(b'\t',3)
                    if len(fields)<3 or int(fields[1])<0 or int(fields[2])<int(fields[1]):
                        raise ValueError('Invalid BED interval')
                    seen += 1
                    if fields[0] in ranks:
                        lines.append(line+b'\n')
                        kept += 1
                    if len(lines)>=chunk_rows:
                        spill()
            if limit and input_rows>=limit:
                break
        if lines:
            spill()
        streams = [p.open('rb') for p in chunks]
        try:
            with output.open('wb') as f:
                f.writelines(heapq.merge(*streams,key=key))
        finally:
            for stream in streams:
                stream.close()
    return seen,kept


def run(args):
    started = time.perf_counter()
    if args.method == 'external':
        seen,kept = external_sort(args.inputs,args.selection,args.output,args.chunk_rows,args.limit)
        emit_metrics(args,dict(project=1,method=args.method,input_rows=seen,output_rows=kept),started)
        return
    _,names,lens = selection(args.selection)
    t = time.perf_counter()
    if args.method in ('radix-threaded','radix-jit-threaded','radix-parallel') and not args.limit and len(args.inputs)>1:
        from concurrent.futures import ThreadPoolExecutor
        with ThreadPoolExecutor(max_workers=min(2,len(args.inputs))) as pool:
            blocks=list(pool.map(lambda path:read_bytes([path]),args.inputs))
        data=np.concatenate(blocks)
        del blocks
    else:
        data = read_bytes(args.inputs,args.limit)
    read_s = time.perf_counter()-t
    t = time.perf_counter()
    records,seen = bed_records(data,names,lens)
    parse_s = time.perf_counter()-t
    t = time.perf_counter()
    order = sorted_order(records,args.method)
    sort_s = time.perf_counter()-t
    t = time.perf_counter()
    if args.method=='radix-parallel':
        from parallel_writer import write_parallel
        write_parallel(args.output,data,records[:,3],records[:,4],order)
    else:
        write_order(args.output,data,records[:,3],records[:,4],order)
    write_s = time.perf_counter()-t
    emit_metrics(args,dict(project=1,method=args.method,input_rows=seen,output_rows=len(order),read_seconds=read_s,parse_seconds=parse_s,sort_seconds=sort_s,write_seconds=write_s),started)


def parser():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--inputs',nargs='+',default=[str(ROOT/'data/shuf.a.bed.gz'),str(ROOT/'data/shuf.b.bed.gz')])
    p.add_argument('--selection',default=str(ROOT/'data/selection.tsv'))
    p.add_argument('--method',choices=METHODS,default='radix-parallel')
    p.add_argument('--sample',default='X')
    p.add_argument('--selection-name',default='provided')
    p.add_argument('--output')
    p.add_argument('--limit',type=int,default=0,help='Input prefix only; 0 processes all records')
    p.add_argument('--chunk-rows',type=int,default=250000)
    p.add_argument('--metrics')
    return p


if __name__ == '__main__':
    args=parser().parse_args()
    if args.limit < 0 or args.chunk_rows < 1:
        raise SystemExit('limit must be nonnegative and chunk-rows positive')
    if args.output is None:
        if any(c in args.sample+args.selection_name for c in '/\\'):
            raise SystemExit('Sample and selection name must be filename components')
        args.output=str(ROOT/'results/sorted_bed_file_per_sample'/f'{args.sample}.{args.selection_name}.sorted.bed')
    run(args)
