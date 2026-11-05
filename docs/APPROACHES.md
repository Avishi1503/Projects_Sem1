# How the four projects work

All implementations are Python. NumPy, pandas and Numba are allowed dependencies in this package; no handwritten C/C++ program or executable compiler is required. If the instructor restricts packages, confirm that separately: "Python only" does not tell us whether JIT compilation is permitted. The ordinary Python and pandas alternatives remain available.

## 1. Quick BED sort

Imagine labelled trays: one tray for each chromosome, arranged exactly as the selection file says. First discard rows whose chromosomes have no tray. Then sort retained rows by **selection rank, numerical start, numerical end**. Equal keys keep their original input order. The statement explicitly specifies chromosome order; sorting coordinates within each chromosome is an additional conventional BED ordering choice, documented here.

Extra BED fields are payload, not sorting keys. One supplied run puts a length in column 4; the other puts a read identifier there. Both formats survive unchanged. Duplicate records are retained, because merging runs is not deduplication.

Eleven alternatives are provided. In-memory sorting comparisons share a compiled parser and serial output writer; `external` uses an independent Python parser and merge path. The explicitly labelled `radix-parallel` additionally tests a parallel output writer.

| Method | Idea | Time after parsing | Main tradeoff |
|---|---|---|---|
| `timsort` | Python sorts a list of row indices using tuple keys | O(N log N) worst case | Easiest explanation; Python objects use substantial memory |
| `bucket` | Partition by chromosome, then Python-sort each bucket | O(N + sum n_c log n_c) | Fewer comparisons; still Python objects |
| `quicksort` | Pack the keys into uint64, use NumPy comparison sort, restore ties | O(N log N) with installed NumPy introsort | Fast general-purpose comparison sort |
| `heapsort` | Same packed keys, NumPy heap sort, restore ties | O(N log N) | Predictable complexity, poorer locality |
| `lexsort` | Stable multi-key numerical sorting | Implementation dependent; integer keys use specialized sorts | Avoids bit packing; handles wide coordinates |
| `radix` | Pack chromosome/start/end into uint64 and use stable integer radix sort | O(wN), fixed-width w | Fast, compact indices; needs temporary arrays |
| `radix-threaded` | Same radix sort, decompress two run files concurrently | Same sorting complexity | Can reduce reading time; extra transient memory and two reader threads |
| `radix-jit` | Explicit stable 8-bit LSD radix passes, written in Python with Numba | O(8N + 8*256) | Moves keys alongside indices for sequential reads; larger temporary arrays |
| `radix-jit-threaded` | Explicit radix passes plus two concurrent file readers | Same sorting complexity | Combines the two optimizations; measured rather than assumed faster |
| `radix-parallel` | Explicit radix, two readers and up to four output-copy workers | Same sorting complexity | Tests whether parallel copying helps the remaining output bottleneck |
| `external` | Sort bounded chunks on disk and heap-merge their streams | O(N log M + N log K) | Low RAM; extra disk writes and Python merge work |

Here M is the chunk size and K the number of chunks. The external version opens K streams; for thousands of chunks a production system should use a multi-stage merge. The supplied data require far fewer.

Packing example: allocate enough bits for every observed start and end, and put chromosome rank in the most significant bits. Comparing one resulting integer is equivalent to comparing the three keys. **Never silently truncate:** if all keys do not fit into 64 bits, this code falls back to `lexsort`. It does not allocate one cell per genomic coordinate; a direct counting sort over billions of positions would waste memory. The explicit radix implementation instead uses just 256 counting buckets per pass. It sorts by the least significant byte first and maintains stability at every pass, so more significant bytes establish the final lexicographic order without breaking earlier ordering.

The built-in quicksort/heapsort names refer to NumPy **1.23.5**, the tested version. Newer NumPy releases may map named options differently. Integer `stable` and `mergesort` are aliases of a specialized radix implementation here, so presenting them as two independent algorithms would be misleading. See the [NumPy sorting documentation](https://numpy.org/doc/stable/reference/generated/numpy.sort.html) and [Python sorting guide](https://docs.python.org/3.13/howto/sorting.html).

Output: `results/sorted_bed_file_per_sample/X.provided.sorted.bed`. Both the sample and selection name appear in the filename. Run `--sample` and `--selection-name` to change them.

## 2. Rescale a fragment-length distribution

Use `end - start` for fragment length. Build a histogram `q[L]` of query rows and normalize the reference weights to `p[L]`. Decide how many rows to retain for each length **before** sampling. Sample original row positions without replacement, then write their complete lines in original order.

If two distinct input records have identical text, both can legitimately appear in the result. "Without replacement" applies to original row positions; it does not silently deduplicate biological reads.

Two explicit policies are included because "matches the reference" and "retains a large dataset" can conflict:

* `capped` (default): request a scale of 30% of input row count, then set `quota[L] = min(q[L], floor(scale * p[L]))`. This yields the intended shape where enough rows exist. Scarce bins are capped and the output can contain fewer than 30% of rows. The 30% is a chosen, configurable setting, not a requirement in the PDF.
* `strict`: on lengths present in both inputs, choose `scale = min(q[L] / p[L])`, then use the same floor quotas. This is the largest continuous proportional scale on shared support, before integer rounding. One scarce length can make the retained sample small. Missing reference bins remain impossible; their mass is explicitly reported. It is not a claim of an exact real-valued probability match after rounding.

No sampling algorithm can create a fragment length absent from the query without fabricating data. Report missing reference mass, retained count, and total variation distance `0.5 * sum(abs(observed_probability-reference_probability))`. A small distance means a closer distribution. Do not choose a "fastest" algorithm just because it quietly writes fewer rows: every method in a policy comparison receives identical quotas.

For this assignment's goal of matching the reference as closely as possible, prefer the **strict output**. The capped output is also provided to show the retention/accuracy tradeoff and uses larger identical quotas for the runtime comparison. The runner creates both, and the report gives their actual counts and distances.

Four sampling approaches:

1. **Sequential exact sampling (Python):** keep the next row with probability `still_needed / still_remaining` within its length bin. Decrement the remaining count each time. This gives an unbiased subset of the requested size.
2. **Sequential exact sampling (Numba):** the same mathematical rule executed as a compiled Python loop. O(N + L), with a compact output index array.
3. **Random priorities:** give each row a random number, group by length and sort by priority, then retain the first quota rows. Uniform subsets, but adds an O(N log N) sorting step and temporary arrays.
4. **Stratified reservoir sampling:** maintain a reservoir of size quota for each length. A later row replaces an existing entry with the appropriate probability. Selection is O(N), but per-row Python RNG calls cost time.

The current implementation first reads and indexes the input for every alternative, so it can validate quotas and preserve rows consistently. The reservoir selector itself only needs O(sum quotas) storage; the delivered pipeline is **not** claimed to be a fully streaming reservoir pipeline. Different generators can select different records with the same seed; all meet identical per-length quotas. Repeating the same method and seed is deterministic in the tested environment.

## 3. V-plot

The provided file has 13 columns: a leading field, a six-column region, and a six-column fragment. One row is an already computed region/fragment overlap. A fragment overlapping two regions contributes once to each region. Do not deduplicate these pairs.

Coordinates used (one-based column numbers):

* Region: chromosome 2, start 3, end 4, strand 7.
* Fragment: chromosome 8, start 9, end 10.
* Fragment length: column 10 minus column 9. The redundant length field is not needed.

Compute the difference between fragment and region midpoints. Mirror it when the **region** strand is negative. To avoid floating-point arithmetic, calculate `doubled_offset = fragment_start + fragment_end - region_start - region_end`, change its sign if needed, then take floor division by 2. Thus a half-base offset of +0.5 falls in bin 0 and -0.5 in bin -1. `--genomic` disables strand mirroring.

The statement does not fully specify midpoint rounding, strand orientation, or the exact last offset column. These are explicit conventions, not facts inferred from a filename. Defaults are offsets **-500 through +499 inclusive** (1,000 bins) and lengths 0 through 1,000. They are configurable. Rows outside the selected window, overly long fragments, and mismatched chromosomes are counted and reported as excluded. The first matrix column contains fragment length; subsequent columns are the offset bins. The long form contains `offset, length, count` and is directly suitable for Gnuplot.

Compare a Python `Counter`, pandas chunks plus `numpy.bincount`, and a Numba byte parser that increments a dense matrix. All are O(N) in row count for fixed matrix size. The differences are parsing overhead and Python object creation. The compiled implementation reads bounded chunks from standard input and does not need the whole file in memory.

## 4. NKX2-1 in TCGA-LUAD

Match every expression file to the sample sheet using **File ID**, check the filename, and include exactly `TCGA-LUAD` / `Primary Tumor` or `Solid Tissue Normal`. Ignore mutation files, other data types and other sample categories. The original ZIP and inner tar.gz are streamed; no multi-gigabyte expression folder must be extracted.

Use the column named `tpm_unstranded`, and match `NKX2-1` exactly in the column named `gene_name`. Transform each TPM as `log2(TPM + 1)`, not `log2(TPM) + 1`. Reject missing/ambiguous gene rows and negative/nonfinite TPM. Check file sizes and MD5 against the supplied manifest. The problem's normal-sample example must produce TPM **53.1718**.

Three parsers are benchmarked with the same archive reading and integrity checks: full pandas table parsing, line-by-line field splitting, and exact byte-token search followed by parsing only candidate rows. Byte search reduces work because only one gene is needed. Archive decompression still has to traverse every file. A cached results table is useful for plotting again, but is **not** counted as a fresh-data benchmark win.

The expression table retains all eligible files. Plots explicitly identify their unit (files or biological samples). For inference, repeated files from the same patient/tissue are aggregated and matched tumor/normal patients are tested with paired Wilcoxon. Unmatched patients can be compared separately with Mann-Whitney. These tests address expression differences, not diagnostic sensitivity, specificity or clinical utility. See the [SciPy Wilcoxon documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.wilcoxon.html) and [Mann-Whitney documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.mannwhitneyu.html).

## How to explain the optimization in a viva

"I established one output definition, wrote edge-case checks, and compared several Python implementations. I measured complete fresh-process runs, including decompression, parsing and writing. I repeated each full-data comparison, reported the median and range, and checked output equivalence. My choices reduce Python work and unnecessary sorting. I report JIT startup separately and do not claim my laptop's winner is universally fastest."

Know these distinctions: algorithm versus implementation; asymptotic complexity versus measured runtime; execution time versus first-use compilation; same random quotas versus identical sampled rows; full data versus screening subsets; descriptive group differences versus paired statistical inference.
