# Benchmark results

The tables below compare the Python implementations on the supplied datasets. Times vary by machine; these are the measurements from the test laptop.

## Fastest methods in these runs

| Project | Fastest tested method | Median execution | Observed range | Peak RSS |
|---|---|---:|---:|---:|
| 1. BED sorting | `radix-parallel` | 14.01 s | 13.21-14.34 s | 1752 MiB |
| 2. Reference resampling | `compiled` | 9.45 s | 9.20-9.60 s | 1222 MiB |
| 3. V-plot matrix | `compiled` | 5.74 s | 5.11-6.10 s | 115 MiB |
| 4. TCGA gene extraction | `bytes` | 15.89 s | 15.80-16.14 s | 69 MiB |

## How the timings were measured

Each run starts a new Python process. The timer includes imports, loading the Numba cache, reading, decompression, calculation and writing the output. Each method uses the same input and output requirements. Runs are repeated at least three times, in a shuffled order with a fixed seed within each batch. The final sorting batch compares the leading methods together. Writes are buffered by the operating system; no forced disk sync is used. Output hashing by the benchmark script happens after timing, while checks inside each program remain included.

The timed outputs are the sorted BED, sampled BED and histogram, V-plot matrix, and TCGA expression table. Plotting and patient-level statistics are measured separately in `benchmarks/figure_times.json`. In particular, the TCGA time in the table does not include making the boxplot.

Disk caches were not cleared and other desktop applications were open. Downloading and installing packages are excluded. Threaded sorting variants use two readers; radix-parallel also uses four output-copy workers. Other variants use one computational worker.

Screening uses smaller inputs: 100,000 rows from each sorting file, 200,000 rows for resampling and V-plot, and 12 eligible TCGA files. Timsort and chromosome buckets were only run at this size because of their Python-object memory use. Full-input results are labelled separately.

## 1. BED sorting

| Scope | Method | Runs | Median seconds | Min-max seconds | Maximum sampled RSS (MiB) |
|---|---|---:|---:|---:|---:|
| full | radix-parallel | 3 | 14.013 | 13.214-14.339 | 1751.8 |
| full | radix-jit-threaded | 6 | 15.693 | 15.355-16.750 | 1751.5 |
| full | radix-jit | 3 | 16.766 | 16.554-17.095 | 1751.4 |
| full | quicksort | 6 | 19.021 | 18.742-21.025 | 1520.5 |
| full | radix-threaded | 3 | 22.677 | 22.025-24.281 | 1513.2 |
| full | radix | 3 | 25.394 | 23.225-25.670 | 1503.3 |
| full | heapsort | 3 | 28.839 | 27.033-29.628 | 1541.8 |
| full | lexsort | 3 | 41.263 | 38.198-42.987 | 1493.8 |
| full | external | 3 | 91.388 | 88.168-93.981 | 131.9 |
| screening | radix-jit-threaded | 3 | 1.741 | 1.694-1.811 | 110.6 |
| screening | radix | 3 | 1.743 | 1.668-1.872 | 111.1 |
| screening | radix-jit | 3 | 1.750 | 1.686-1.817 | 111.6 |
| screening | lexsort | 3 | 1.792 | 1.748-1.875 | 111.0 |
| screening | quicksort | 3 | 1.805 | 1.695-1.816 | 111.6 |
| screening | heapsort | 3 | 1.814 | 1.748-1.865 | 108.3 |
| screening | radix-threaded | 3 | 1.816 | 1.758-2.919 | 110.1 |
| screening | radix-parallel | 3 | 1.876 | 1.809-1.884 | 109.6 |
| screening | external | 3 | 1.922 | 1.871-1.942 | 111.5 |
| screening | bucket | 3 | 1.996 | 1.932-2.057 | 114.6 |
| screening | timsort | 3 | 2.109 | 1.989-2.134 | 136.8 |

Input rows: **13,835,737**. Retained after chromosome selection: **13,452,140**. All full-data outputs equal: **True**.

## 2. Reference resampling

| Scope | Method | Runs | Median seconds | Min-max seconds | Maximum sampled RSS (MiB) |
|---|---|---:|---:|---:|---:|
| full | compiled | 3 | 9.448 | 9.195-9.603 | 1221.6 |
| full | priority | 3 | 29.576 | 29.363-30.238 | 1410.4 |
| full | python | 3 | 50.874 | 45.932-53.684 | 1329.6 |
| full | reservoir | 3 | 52.375 | 50.855-52.746 | 1217.4 |
| screening | priority | 3 | 2.007 | 1.802-2.128 | 99.2 |
| screening | reservoir | 3 | 2.425 | 2.423-2.434 | 97.1 |
| screening | python | 3 | 2.489 | 2.372-2.564 | 98.9 |
| screening | compiled | 3 | 2.873 | 1.799-3.067 | 110.1 |

Query rows: **13,835,737**. Capped-policy sampled rows: **4,026,230**. Missing reference mass: **0.000000%**. Total variation from reference: **0.027844**. Every selector meets the same quotas and uses distinct input indices.

Random sample records can differ between methods; their per-length counts must agree. `results/rescale/query.strict.summary.json` documents the stricter shape/retention tradeoff separately.

## 3. V-plot matrix

| Scope | Method | Runs | Median seconds | Min-max seconds | Maximum sampled RSS (MiB) |
|---|---|---:|---:|---:|---:|
| full | compiled | 3 | 5.739 | 5.110-6.104 | 115.0 |
| full | pandas | 3 | 8.142 | 7.765-8.181 | 141.5 |
| full | dictionary | 3 | 19.037 | 18.154-19.500 | 104.3 |
| screening | dictionary | 3 | 2.265 | 2.251-2.364 | 82.8 |
| screening | compiled | 3 | 2.364 | 2.361-2.387 | 126.9 |
| screening | pandas | 3 | 2.377 | 2.316-2.435 | 136.5 |

Input region/fragment pairs: **4,336,134**. Counted in the configured window: **4,234,419**. Excluded: **101,715**. All full-data matrices equal: **True**.

## 4. TCGA gene extraction

| Scope | Method | Runs | Median seconds | Min-max seconds | Maximum sampled RSS (MiB) |
|---|---|---:|---:|---:|---:|
| full | bytes | 3 | 15.893 | 15.797-16.143 | 68.9 |
| full | lines | 3 | 37.064 | 36.004-37.105 | 69.1 |
| full | pandas | 3 | 54.050 | 52.866-54.636 | 110.7 |
| screening | bytes | 3 | 1.125 | 1.119-1.190 | 63.3 |
| screening | lines | 3 | 1.630 | 1.433-1.630 | 66.2 |
| screening | pandas | 3 | 2.627 | 2.555-2.689 | 105.0 |

Eligible expression files: **598**. Manifest MD5 checks: **598**. Groups: `{'Primary Tumor': 539, 'Solid Tissue Normal': 59}`. All extraction tables equal: **True**.

## First-use compilation cost

The main timings use an existing Numba cache. This smaller experiment compares the first call with an empty cache, the first call with a saved cache, and a second call in the same process. The first-call times include a small example calculation as well as compilation or loading.

| Kernel | Empty-cache first call | Populated-cache first call | In-process repeat |
|---|---:|---:|---:|
| BED parser | 2.1742 s | 0.5443 s | 0.000004 s |
| length parser | 0.8191 s | 0.0321 s | 0.000004 s |
| tie ordering | 1.3594 s | 0.0323 s | 0.000015 s |
| exact sampling | 0.6625 s | 0.0350 s | 0.000017 s |
| V-plot | 1.4003 s | 0.0340 s | 0.001608 s |
| BED output copy | 0.3814 s | 0.0367 s | 0.000017 s |
| Parallel BED copy | 0.7914 s | 0.0668 s | 0.000019 s |
| Explicit radix sort | 1.3058 s | 0.0443 s | 0.000016 s |

Run `python benchmarks/jit_cost.py` to repeat this check. A first run on another computer includes compilation. NumPy and pandas use libraries compiled when installed.

## Complete first-use runs

These full-input runs use a separate, empty Numba cache and include compilation. Each time is a single observation, not a median. The operating-system file cache was still uncontrolled. These runs are kept separate from the repeated timings above.

| Project | Selected method | Empty-cache complete run |
|---|---|---:|
| 1 | radix-parallel | 15.71 s |
| 2 | compiled | 10.66 s |
| 3 | compiled | 7.23 s |

TCGA byte extraction uses no Numba kernel; its main fresh-process timings already include Python/library startup.

## Recommended resampling output

For the closest reference match, use `results/rescale/query.strict.bed`: **1,220,621 rows**, total variation **0.00018906**. The capped comparison output retains more rows but has greater shape error. The two policies produce different amounts of data, so their runtimes are not directly comparable.

## Expression interpretation

The files contain 598 eligible expression records, 587 biological samples and 517 patients. The descriptive boxplot shows each expression file. Patient-level testing uses 58 matched tumor/normal pairs after aggregating repeated observations.
Paired Wilcoxon p-value: **0.00383082**. Median paired change in log2(TPM + 1): **-0.3113**.

The paired test compares tumor and normal tissue within patients. The all-file boxplot shows the overall spread. Neither result alone establishes diagnostic accuracy.

Across all files, the tumor and normal medians are close, but tumor values vary more. Among matched patients, the median tumor-minus-normal change is negative. Only one normal patient is unmatched, so a separate unmatched-group test is not reported. The example plot in the assignment is illustrative; its p-value is not a target for this analysis.

## Supporting files

* [Algorithm notes](docs/APPROACHES.md)
* [Setup and run commands](README.md)
* Raw observations: `benchmarks/raw_results.jsonl`
* Full command logs: `benchmarks/logs/`
* Per-stage timings: `benchmarks/summary.csv`
* Input provenance and SHA-256: `data/sources.json`
* Source hashes from the original timing runs: `benchmarks/source_hashes.json` (before the formatting and comment edits)

Tested environment: Python 3.10.9, NumPy 1.23.5, pandas 1.5.3, Numba 0.56.4, SciPy 1.10.0; Intel Core i5-1135G7, 4 physical cores / 8 logical cores, 15.7 GiB RAM.
