"""Project 3: mapped BED on stdin -> dense or long V-plot count matrix."""
import argparse
import collections
import io
import time
from pathlib import Path
import numpy as np
from common import ROOT,open_input,emit_metrics,json_write
from kernels import vplot_rows

METHODS=['compiled','dictionary','pandas']


def calculate(stream,method,low,high,max_length,oriented=True,limit=0):
    matrix=np.zeros((max_length+1,high-low+1),np.int64)
    seen=kept=mismatched=0
    if method=='compiled':
        while True:
            if limit:
                lines=[]
                for _ in range(min(50000,limit-seen)):
                    line=stream.readline()
                    if not line: break
                    lines.append(line)
                raw=b''.join(lines)
            else:
                raw=stream.read(8*1024*1024)
                if raw and not raw.endswith(b'\n'):
                    raw+=stream.readline()
            if not raw: break
            part,n,k,mis=vplot_rows(np.frombuffer(raw,np.uint8),low,high,max_length,oriented)
            matrix+=part; seen+=n; kept+=k; mismatched+=mis
            if limit and seen>=limit: break
    elif method=='dictionary':
        counts=collections.Counter()
        for line in stream:
            if not line.strip() or line.startswith(b'#'): continue
            f=line.rstrip(b'\r\n').split(b'\t')
            if len(f)!=13: raise ValueError('Expected 13 mapped BED columns')
            seen+=1
            rs,re,fs,fe=map(int,(f[2],f[3],f[8],f[9]))
            if min(rs,re,fs,fe)<0 or re<rs or fe<fs: raise ValueError('Invalid interval')
            if f[1]!=f[7]: mismatched+=1
            else:
                d=fs+fe-rs-re
                if oriented and f[6]==b'-': d=-d
                offset=d//2; length=fe-fs
                if low<=offset<=high and length<=max_length:
                    counts[length,offset-low]+=1; kept+=1
            if limit and seen>=limit: break
        for (length,column),count in counts.items(): matrix[length,column]=count
    else:
        import pandas as pd
        for frame in pd.read_csv(stream,sep='\t',header=None,comment='#',usecols=[1,2,3,6,7,8,9],chunksize=100000,nrows=limit or None):
            rs=frame[2].to_numpy(dtype=np.int64); re=frame[3].to_numpy(dtype=np.int64)
            fs=frame[8].to_numpy(dtype=np.int64); fe=frame[9].to_numpy(dtype=np.int64)
            if np.any(rs<0)|np.any(fs<0)|np.any(re<rs)|np.any(fe<fs): raise ValueError('Invalid interval')
            length=fe-fs; delta=fs+fe-rs-re
            if oriented: delta=np.where(frame[6].to_numpy()=='-',-delta,delta)
            offset=delta//2
            same=frame[1].to_numpy()==frame[7].to_numpy()
            mask=same&(offset>=low)&(offset<=high)&(length<=max_length)
            flat=length[mask]*(high-low+1)+offset[mask]-low
            matrix+=np.bincount(flat,minlength=matrix.size).reshape(matrix.shape)
            seen+=len(frame); kept+=int(mask.sum()); mismatched+=int((~same).sum())
    return matrix,dict(input_rows=seen,counted_rows=kept,excluded_rows=seen-kept,chromosome_mismatches=mismatched,offset_low=low,offset_high=high,max_length=max_length,orientation='region-strand' if oriented else 'genomic',rounding='floor signed doubled-center difference / 2')


def write_matrix(path,matrix,low,long_form=False):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    if long_form:
        with path.open('w',encoding='ascii',newline='\n') as f:
            for length,row in enumerate(matrix):
                for j,count in enumerate(row): f.write(f'{low+j}\t{length}\t{count}\n')
                f.write('\n')
    else:
        np.savetxt(path,np.column_stack((np.arange(len(matrix)),matrix)),fmt='%d',delimiter='\t')
    np.save(path.with_suffix('.npy'),matrix)


def run(args):
    started=time.perf_counter()
    with open_input(args.input) as f:
        matrix,info=calculate(f,args.method,args.low,args.high,args.max_length,not args.genomic,args.limit)
    calc_s=time.perf_counter()-started
    t=time.perf_counter(); write_matrix(args.output,matrix,args.low,args.long); write_s=time.perf_counter()-t
    if int(matrix.sum())!=info['counted_rows']: raise AssertionError('Count conservation failed')
    info.update(project=3,method=args.method,calculate_seconds=calc_s,write_seconds=write_s,matrix_shape=list(matrix.shape),matrix_sum=int(matrix.sum()),format='long' if args.long else 'wide')
    json_write(Path(args.output).with_suffix('.summary.json'),info)
    emit_metrics(args,info,started)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',default='-',help='- reads decompressed BED from stdin; .gz is accepted as a file')
    p.add_argument('--output',default=str(ROOT/'results/vplot/vplot.matrix.tsv'))
    p.add_argument('--method',choices=METHODS,default='compiled')
    p.add_argument('--low',type=int,default=-500)
    p.add_argument('--high',type=int,default=499)
    p.add_argument('--max-length',type=int,default=1000)
    p.add_argument('--genomic',action='store_true',help='Disable reverse-strand mirroring')
    p.add_argument('--long',action='store_true')
    p.add_argument('--limit',type=int,default=0)
    p.add_argument('--metrics')
    args=p.parse_args()
    if args.low>args.high or args.max_length<0 or args.limit<0: p.error('Invalid matrix range or limit')
    if (args.high-args.low+1)*(args.max_length+1)>20_000_000: p.error('Matrix too large; reduce ranges')
    run(args)
