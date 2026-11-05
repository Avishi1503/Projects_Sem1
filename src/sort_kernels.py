"""Explicit stable LSD radix sort, implemented in Python and compiled by Numba."""
import numpy as np
from numba import njit


@njit(cache=True)
def radix_indices(keys):
    # Move keys alongside row indices, so every subsequent key scan is linear.
    # This avoids randomly gathering the original keys on every radix pass.
    n=len(keys)
    current_keys=keys.copy()
    next_keys=np.empty(n,np.uint64)
    current_indices=np.arange(n,dtype=np.int64)
    next_indices=np.empty(n,np.int64)
    counts=np.zeros(256,np.int64)
    for shift in range(0,64,8):
        counts[:]=0
        for i in range(n): counts[int((current_keys[i]>>np.uint64(shift))&np.uint64(255))]+=1
        if counts.max()==n: continue
        total=0
        for bucket in range(256):
            count=counts[bucket]; counts[bucket]=total; total+=count
        for i in range(n):
            key=current_keys[i]
            bucket=int((key>>np.uint64(shift))&np.uint64(255))
            dest=counts[bucket]; counts[bucket]+=1
            next_keys[dest]=key
            next_indices[dest]=current_indices[i]
        current_keys,next_keys=next_keys,current_keys
        current_indices,next_indices=next_indices,current_indices
    return current_indices
