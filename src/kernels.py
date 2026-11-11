"""Numba routines for reading BED rows, sampling and counting."""
import numpy as np
from numba import njit


@njit(cache=True)
def has_prefix(b, begin, stop, prefix):
    if stop - begin < len(prefix):
        return False
    for j in range(len(prefix)):
        if b[begin + j] != prefix[j]:
            return False
    return True


@njit(cache=True)
def uint_field(b, start, stop):
    if stop <= start:
        raise ValueError("Empty integer field")
    x = np.uint64(0)
    for k in range(start, stop):
        c = b[k]
        if c < 48 or c > 57:
            raise ValueError("Expected a nonnegative integer")
        if x > np.uint64(1844674407370955161) or (x == np.uint64(1844674407370955161) and c > 53):
            raise ValueError("Coordinate exceeds uint64")
        x = x * np.uint64(10) + np.uint64(c - 48)
    return x


@njit(cache=True)
def bed_records(b, names, name_lengths):
    """Columns: selection rank, start, end, byte offset, complete line length."""
    n = 0
    for c in b:
        if c == 10:
            n += 1
    if len(b) and b[-1] != 10:
        n += 1
    result = np.empty((n, 5), np.uint64)
    pos = 0
    used = 0
    seen = 0
    while pos < len(b):
        begin = pos
        while pos < len(b) and b[pos] != 10:
            pos += 1
        stop = pos
        if stop > begin and b[stop - 1] == 13:
            stop -= 1
        line_end = min(pos + 1, len(b))
        pos = line_end
        if stop == begin or b[begin] == 35:
            continue
        # Skip optional BED headers.
        if has_prefix(b, begin, stop, b"track "):
            continue
        if has_prefix(b, begin, stop, b"browser "):
            continue
        t1 = begin
        while t1 < stop and b[t1] != 9:
            t1 += 1
        t2 = t1 + 1
        while t2 < stop and b[t2] != 9:
            t2 += 1
        t3 = t2 + 1
        while t3 < stop and b[t3] != 9:
            t3 += 1
        if t1 == stop or t2 >= stop:
            raise ValueError("BED requires three tab-separated columns")
        start = uint_field(b, t1 + 1, t2)
        end = uint_field(b, t2 + 1, t3)
        if end < start:
            raise ValueError("BED end precedes start")
        seen += 1
        rank = -1
        for j in range(len(name_lengths)):
            if name_lengths[j] == t1 - begin:
                equal = True
                for k in range(name_lengths[j]):
                    if b[begin + k] != names[j, k]:
                        equal = False
                        break
                if equal:
                    rank = j
                    break
        if rank >= 0:
            result[used, 0] = rank
            result[used, 1] = start
            result[used, 2] = end
            result[used, 3] = begin
            result[used, 4] = line_end - begin
            used += 1
    return result[:used], seen


@njit(cache=True)
def bed_offsets_lengths(b):
    n = 0
    for c in b:
        n += c == 10
    if len(b) and b[-1] != 10:
        n += 1
    result = np.empty((n, 3), np.int64)
    pos = 0
    used = 0
    while pos < len(b):
        begin = pos
        while pos < len(b) and b[pos] != 10:
            pos += 1
        stop = pos
        pos = min(pos + 1, len(b))
        if stop > begin and b[stop - 1] == 13:
            stop -= 1
        if stop == begin or b[begin] == 35:
            continue
        if has_prefix(b, begin, stop, b"track "):
            continue
        if has_prefix(b, begin, stop, b"browser "):
            continue
        t1 = begin
        while t1 < stop and b[t1] != 9:
            t1 += 1
        t2 = t1 + 1
        while t2 < stop and b[t2] != 9:
            t2 += 1
        t3 = t2 + 1
        while t3 < stop and b[t3] != 9:
            t3 += 1
        if t1 == stop or t2 >= stop:
            raise ValueError("BED requires three tab-separated columns")
        start = uint_field(b, t1 + 1, t2)
        end = uint_field(b, t2 + 1, t3)
        if end < start or end - start > 2147483647:
            raise ValueError("Invalid or excessively large fragment length")
        result[used, 0] = begin
        result[used, 1] = pos - begin
        result[used, 2] = end - start
        used += 1
    return result[:used]


@njit(cache=True)
def stable_ties(order, keys):
    """Restore input order within equal keys after unstable comparison sorts."""
    i = 0
    while i < len(order):
        j = i + 1
        while j < len(order) and keys[order[j]] == keys[order[i]]:
            j += 1
        if j - i > 1:
            order[i:j].sort()
        i = j
    return order


@njit(cache=True)
def gather(b, offsets, lengths, order, begin, buffer):
    used = 0
    i = begin
    while i < len(order):
        idx = order[i]
        at = int(offsets[idx])
        size = int(lengths[idx])
        extra = 0 if size and b[at + size - 1] == 10 else 1
        if size + extra > len(buffer):
            raise ValueError("BED line exceeds output buffer")
        if used + size + extra > len(buffer):
            break
        for k in range(size):
            buffer[used + k] = b[at + k]
        used += size
        if extra:
            buffer[used] = 10
            used += 1
        i += 1
    return i, used


@njit(cache=True)
def exact_select(lengths, counts, quotas, seed):
    """Uniform subset per length: P(keep next) = needed / remaining."""
    np.random.seed(seed)
    remaining = counts.copy()
    need = quotas.copy()
    selected = np.empty(int(quotas.sum()), np.int64)
    used = 0
    for i in range(len(lengths)):
        k = lengths[i]
        if need[k] > 0:
            if need[k] == remaining[k] or np.random.random() * remaining[k] < need[k]:
                selected[used] = i
                used += 1
                need[k] -= 1
        remaining[k] -= 1
    if used != len(selected):
        raise ValueError("Sampling failed to meet quotas")
    return selected


@njit(cache=True)
def vplot_rows(b, low, high, max_length, oriented):
    """13-column mapped BED: region cols 2-4/7, fragment cols 8-10.

    Offset = floor(strand_sign * (fragment_start+end-region_start-end)/2).
    The input is already an overlap join: each row is one fragment-region pair.
    """
    matrix = np.zeros((max_length + 1, high - low + 1), np.int64)
    fields = np.empty((13, 2), np.int64)
    pos = 0
    seen = kept = mismatched = 0
    while pos < len(b):
        begin = pos
        while pos < len(b) and b[pos] != 10:
            pos += 1
        stop = pos
        pos = min(pos + 1, len(b))
        if stop > begin and b[stop - 1] == 13:
            stop -= 1
        if begin == stop or b[begin] == 35:
            continue
        col = 0
        at = begin
        for k in range(begin, stop + 1):
            if k == stop or b[k] == 9:
                if col < 13:
                    fields[col, 0] = at
                    fields[col, 1] = k
                at = k + 1
                col += 1
        if col != 13:
            raise ValueError("mapped BED must have exactly 13 columns")
        rs = int(uint_field(b, fields[2, 0], fields[2, 1]))
        re = int(uint_field(b, fields[3, 0], fields[3, 1]))
        fs = int(uint_field(b, fields[8, 0], fields[8, 1]))
        fe = int(uint_field(b, fields[9, 0], fields[9, 1]))
        if re < rs or fe < fs:
            raise ValueError("Reversed interval in mapped BED")
        a, ae = fields[1, 0], fields[1, 1]
        c, ce = fields[7, 0], fields[7, 1]
        same = ae - a == ce - c
        if same:
            for j in range(ae - a):
                if b[a + j] != b[c + j]:
                    same = False
                    break
        seen += 1
        if not same:
            mismatched += 1
            continue
        delta2 = fs + fe - rs - re
        if oriented and b[fields[6, 0]] == 45:
            delta2 = -delta2
        offset = delta2 // 2
        length = fe - fs
        if low <= offset <= high and 0 <= length <= max_length:
            matrix[length, offset - low] += 1
            kept += 1
    return matrix, seen, kept, mismatched
