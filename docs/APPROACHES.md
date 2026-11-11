# Algorithm notes

The four tasks spend time on different things. Sorting needs an ordering algorithm; resampling needs a way to choose rows; a V-plot needs counts; gene extraction needs to find one value in each file. Using sorting everywhere would add work to the last three tasks.

All implementations are Python, with NumPy, pandas and Numba where useful. The fastest methods for the first three tasks use Numba compilation. If packages or JIT compilation are restricted for assessment, that changes which methods can be used. Even the ordinary sorting and resampling variants share some Numba parsing and writing routines.

## 1. BED sorting

Read both BED files, discard chromosomes absent from `selection.tsv`, and sort the remaining rows by:

1. Chromosome position in the selection file.
2. Start coordinate.
3. End coordinate.

Rows with equal keys keep their original order. Extra fields and duplicate rows are kept. Chromosome order is specified by the assignment; sorting by coordinates within a chromosome is the convention used here.

### Methods compared

| Method | How it works | Main cost or limitation |
|---|---|---|
| `timsort` | Python `sorted()` with tuple keys | O(N log N) worst case; many Python objects |
| `bucket` | Group by chromosome, then sort each group | O(N + sum n_c log n_c); still uses Python objects |
| `quicksort` | NumPy sorts packed integer keys, then equal-key ties are restored | Comparison sorting plus a tie pass |
| `heapsort` | Same keys, using NumPy heapsort | O(N log N); slower in these measurements |
| `lexsort` | Sort the three numerical columns directly | Avoids packing; also used if keys exceed 64 bits |
| `radix` | NumPy's stable sort on unsigned integer keys | Fixed-width integer sorting |
| `radix-threaded` | Same sort, with two file-reading threads | Tests whether concurrent decompression helps |
| `radix-jit` | Explicit least-significant-byte radix passes with Numba | Temporary key and index arrays |
| `radix-jit-threaded` | Explicit radix sort plus two readers | Extra memory during reading |
| `radix-parallel` | Explicit radix sort, two readers and four output-copy workers | Fastest here, but needs about 1.7 GiB RAM |
| `external` | Sort chunks on disk, then merge them with a heap | Lower RAM use, more disk I/O |

Most in-memory methods share parsing and writing code. `radix-parallel` changes the output writer too; `external` has its own Python parsing and merging path. These timings therefore compare complete implementations, not just sort functions.

### Why radix sort helps

The chromosome rank, start and end are packed into one unsigned 64-bit integer. The chromosome rank occupies the highest bits, followed by start and end. Sorting these integers gives the same order as comparing the three fields. If the values need more than 64 bits, the program falls back to `lexsort`.

The explicit radix version handles one byte at a time, starting with the lowest byte. Each pass counts 256 possible values and moves rows into the corresponding positions. Each pass is stable, so sorting the next byte preserves the order established by earlier bytes. With eight bytes, its cost is O(8N + 8*256).

Keys move alongside their row indices, making the next scan sequential. The parallel writer computes each row's destination before copying, so workers write to separate slices of the output buffer. Only the copy step is parallel; the file is written in order.

For low-memory machines, `external --chunk-rows 250000` is an alternative. With chunk size M and K chunks, sorting and merging cost O(N log M + N log K). The implementation opens K streams; very large K would need a multi-stage merge.

NumPy's sort behavior depends on the version and data type. In the tested NumPy 1.23.5 environment, stable integer sorting uses radix sort. `stable` and `mergesort` are therefore not counted as separate algorithms here.

## 2. Fragment resampling

Fragment length is `end - start`. Let `q[L]` be the number of query rows of length L and `p[L]` the normalized reference probability. Calculate a quota for each length, then sample that many original row positions without replacement.

The two quota policies are:

- **Capped:** use `scale = fraction * total_rows`, then `quota[L] = min(q[L], floor(scale * p[L]))`. The default fraction is 0.30, a configurable choice rather than an assignment requirement. Scarce lengths are capped, so the result may contain fewer than 30% of the rows.
- **Strict:** on lengths present in both inputs, use `scale = min(q[L] / p[L])`, then the same floor quotas. This gives the largest proportional scale on the shared lengths before integer rounding. A scarce length can limit the whole sample.

Use the **strict output** for the closest reference match. It retains 1,220,621 rows with total variation distance about 0.000189. The capped output retains 4,026,230 rows with distance about 0.027844. The distance is `0.5 * sum(abs(sample_probability - reference_probability))`; smaller is closer.

A length missing from the query cannot be sampled. Its reference probability is reported as missing mass. Integer rounding also prevents a perfectly exact probability match in general.

### Selecting the rows

The Python and Numba sequential methods use the same rule: keep the next row with probability `still_needed / still_remaining` for its length. This gives a uniform subset of the required size. Numba runs that loop with less Python overhead, giving O(N + L) selection work.

The other methods are random priorities (assign a random value, sort within each length, take the first quota rows) and reservoir sampling (maintain a fixed-size sample for each length). Priorities add a sorting step; the reservoir version makes many Python-level random calls.

All methods receive identical quotas. Different random generators may choose different rows, even with the same seed, but the length counts must agree. The output follows input order. Two identical input lines can both be selected because they are distinct original records.

All four implementations read and index the input first. Although reservoir selection can use bounded memory, this particular program is not a fully streaming reservoir implementation.

## 3. V-plot

Each row in `mapped.bed.gz` is an existing region/fragment overlap. A fragment overlapping two regions contributes once to each region.

Using one-based column numbers:

| Field | Columns |
|---|---|
| Region chromosome, start, end, strand | 2, 3, 4, 7 |
| Fragment chromosome, start, end | 8, 9, 10 |

Calculate fragment length and the difference between the fragment and region midpoints. To avoid floating-point arithmetic:

```text
doubled_offset = fragment_start + fragment_end - region_start - region_end
if region_strand is negative: doubled_offset = -doubled_offset
offset = doubled_offset // 2
```

Floor division puts +0.5 in bin 0 and -0.5 in bin -1. `--genomic` disables strand mirroring. These conventions are documented because the assignment does not fully specify rounding, orientation or the last offset column.

The default matrix has lengths 0 through 1,000 and offsets -500 through +499. Rows outside the window, over the length limit, or on mismatched chromosomes are counted as excluded. The wide TSV starts with a length column; the long TSV contains `offset, length, count` for Gnuplot.

Three methods were compared: a Python `Counter`, pandas chunks with `numpy.bincount`, and a Numba byte parser that increments the matrix directly. Each is O(N) for fixed matrix size. The Numba method avoids creating Python objects for every field and reads bounded chunks, so it uses little memory. No sorting is needed.

## 4. NKX2-1 expression

The program streams the ZIP and its inner tar.gz archive. Each expression file is matched to the sample sheet by File ID and filename. Only TCGA-LUAD primary tumor and solid tissue normal samples are included. File sizes and MD5 hashes are checked against the manifest.

The target is an exact `NKX2-1` match in `gene_name`, with expression read from the named `tpm_unstranded` column. The plotted value is `log2(TPM + 1)`. Missing or duplicate gene matches and invalid TPM values are rejected. The example file in the assignment gives TPM 53.1718.

The three approaches parse a pandas table, split each line, or search for the exact gene token in the file's bytes. Byte search was fastest because only one gene is needed. It still decompresses the archive and checks candidate rows against the correct column.

The output has 598 expression files: 539 tumor and 59 normal. The all-file boxplot is descriptive. For statistical testing, repeated observations are aggregated by biological sample and then patient/tissue; the 58 matched patients are compared with a paired Wilcoxon test. The resulting p-value is about 0.00383. Only one normal patient remains unmatched, so a separate unmatched-group test is not reported.

The matched patients show a negative median tumor-minus-normal change, even though the two pooled medians are close. These are different comparisons. An expression difference by itself does not establish diagnostic accuracy.

## Points to understand before presenting

- Big-O describes growth with input size; parsing, allocation and disk I/O also affect measured time.
- First-use Numba compilation is separate from later runs using its cache.
- Equal sampling quotas do not mean identical randomly chosen rows.
- The sorting variants use different amounts of memory and different numbers of workers.
- The timings include file reading and writing. Plotting is measured separately.
- The reported winner is the fastest tested on this machine, not a guarantee for every input or computer.

See the [full results](../REPORT.md) for timing ranges and the [README](../README.md#references) for library references.
