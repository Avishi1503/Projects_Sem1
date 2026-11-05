"""Parallel, disjoint byte copies for sorted BED output; Python/Numba only."""
from pathlib import Path
import numpy as np
import numba
from numba import njit,prange


@njit(parallel=True,cache=True)
def copy_rows(data,offsets,lengths,order,prefix,buffer):
    for j in prange(len(order)):
        row=order[j]
        source=int(offsets[row]); size=int(lengths[row]); dest=prefix[j]
        for k in range(size): buffer[dest+k]=data[source+k]


def write_parallel(path,data,offsets,lengths,order,chunk_rows=150000):
    # read_bytes has already appended a newline to each input's final record.
    # Each worker owns disjoint output slices, so no locking or atomic adds.
    numba.set_num_threads(min(4,numba.config.NUMBA_NUM_THREADS))
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('wb') as f:
        for begin in range(0,len(order),chunk_rows):
            rows=order[begin:begin+chunk_rows]
            prefix=np.empty(len(rows)+1,np.int64); prefix[0]=0
            np.cumsum(lengths[rows],dtype=np.int64,out=prefix[1:])
            buffer=np.empty(prefix[-1],np.uint8)
            copy_rows(data,offsets,lengths,rows,prefix,buffer)
            f.write(memoryview(buffer))
