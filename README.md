# Python bioinformatics projects: implementations and measured comparisons

Start with **[START_HERE.html](START_HERE.html)** for the visual results, **[REPORT.md](REPORT.md)** for measured timings, and **[docs/APPROACHES.md](docs/APPROACHES.md)** for the algorithm explanation and viva notes.

The package covers all four supplied tasks. The goal is fast **execution in Python**. It includes alternatives, complete output files, correctness checks and repeatable benchmarks. The fastest method means the fastest **among those actually tested on this computer**, with stated memory and startup costs.

## Run on this computer

Open PowerShell in this folder. The existing Python environment already contains the needed packages:

```powershell
& 'C:\ProgramData\anaconda3\python.exe' run_all.py
```

Run just one project:

```powershell
& 'C:\ProgramData\anaconda3\python.exe' run_all.py --projects 1
& 'C:\ProgramData\anaconda3\python.exe' run_all.py --projects 2 3 4
```

The runner reads the measured winners from `benchmarks/winners.json`. It overwrites the generated results for the selected projects. The full data and original supplied documents are in this folder; no download is needed. The original files in Downloads remain untouched.

## Run on another computer

Use Python 3.10 and install the recorded dependencies into an isolated environment:

```text
python -m venv .venv
# Activate .venv for your operating system, then:
python -m pip install -r requirements-tested.txt
python run_all.py
```

The exact package versions record the tested environment. These older pins target Python 3.10, not Python 3.13. There is no C/C++ source to compile. Numba compiles selected Python loops on first use and saves a local cache; first-run execution can therefore be slower. Source code and numerical arrays are still authored entirely in Python.

If you received the small submission ZIP, copy the supplied archives to a local directory and prepare them:

```text
python prepare_data.py --downloads /path/to/local/downloads
```

This reads `27087241.zip`, `mapped.bed.gz` and `27129882.zip`, prepares their required local contents and records SHA-256 hashes. It has no download code. Input locations are also configurable on each project script. The runner's `--input-dir` accepts a complete prepared data folder.

## Commands for each project

From the project folder, `python` below means your chosen Python 3.10 environment. On this PC replace it with `& 'C:\ProgramData\anaconda3\python.exe'` in PowerShell.

### Project 1

```text
python src/sort_bed.py --inputs data/shuf.a.bed.gz data/shuf.b.bed.gz --selection data/selection.tsv --sample X --selection-name provided --method radix-parallel
```

Alternatives: `radix`, `radix-threaded`, `radix-jit`, `radix-jit-threaded`, `radix-parallel`, `quicksort`, `heapsort`, `lexsort`, `timsort`, `bucket`, `external`. The default output is `results/sorted_bed_file_per_sample/X.provided.sorted.bed`. Chromosome order comes from the selection file; ties are resolved by start, end, then original input order. All extra fields and duplicates are retained. `external --chunk-rows 250000` is the limited-RAM alternative. Timsort and bucket methods can consume several GB of Python object memory on the complete input; use `--limit` for experimentation on small prefixes.

### Project 2

```text
python src/rescale.py --input data/query.bed.gz --reference data/reference.hist --method compiled --policy capped --fraction 0.30 --seed 2026
python src/rescale.py --input data/query.bed.gz --reference data/reference.hist --method compiled --policy strict --output results/rescale/query.strict.bed
python src/make_figures.py --project 2
```

Alternatives: `python`, `priority`, `reservoir`. Every method samples without replacement and exactly meets the same integer quotas. **Use `query.strict.bed` when the priority is the closest reference match.** The capped policy trades shape for retention when bins are scarce, and its larger quotas are used for the algorithm runtime comparison. The strict policy uses the largest continuous proportional scale on shared support, then rounds down. The runner creates both outputs. Compare their JSON summaries and histograms: this is a scientific-output choice, separate from execution speed.

### Project 3

```text
python src/vplot.py --input data/mapped.bed.gz --method compiled
python src/vplot.py --input data/mapped.bed.gz --method compiled --long --output results/vplot/vplot.long.tsv
python src/make_figures.py --project 3
```

Alternatives: `dictionary`, `pandas`. Wide output has a fragment-length column followed by 1,000 offset columns for -500 through +499. The default length range is 0 through 1,000. Configure with `--low`, `--high`, `--max-length`; `--genomic` disables region-strand mirroring. Exclusions are reported rather than silently lost.

Standard input is supported as requested. On Linux/macOS:

```sh
gzip -dc data/mapped.bed.gz | python src/vplot.py --input - --method compiled
```

On Windows use `--input data/mapped.bed.gz` to avoid PowerShell text-pipe encoding conversions. For a portable stdin example, this Python wrapper streams decompressed bytes:

```python
import gzip, shutil, subprocess, sys
p = subprocess.Popen([sys.executable, 'src/vplot.py', '--input', '-'], stdin=subprocess.PIPE)
with gzip.open('data/mapped.bed.gz', 'rb') as f:
    shutil.copyfileobj(f, p.stdin)
p.stdin.close()
assert p.wait() == 0
```

The Gnuplot recipe is `docs/vplot.gnuplot`, used after generating the long matrix. The delivered PNG is rendered with Matplotlib; the Gnuplot executable was not available for a local rendering test.

### Project 4

```text
python src/tcga.py --archive data/tcga_source.zip --gene NKX2-1 --method bytes
python src/make_figures.py --project 4
```

Alternatives: `lines`, `pandas`. All read the same archive and verify each eligible expression file against the manifest. The optimized method finds the exact gene token, then reads the named TPM field. Change `--gene` and `--output` to extract another gene; the supplied figure script specifically visualizes NKX2-1. The archive is streamed, saving the space and time of fully extracting it.

## Reproduce the comparisons

```text
python tests/test_correctness.py
python benchmarks/benchmark.py --project 1 --scope screening --repeats 3
python benchmarks/benchmark.py --project 1 --scope full --methods radix radix-threaded radix-jit radix-jit-threaded radix-parallel quicksort heapsort lexsort external --repeats 3
python benchmarks/benchmark.py --project 2 --scope full --repeats 3
python benchmarks/benchmark.py --project 3 --scope full --repeats 3
python benchmarks/benchmark.py --project 4 --scope full --repeats 3
python benchmarks/jit_cost.py
python benchmarks/cold_start_full.py
python benchmarks/report.py
```

The harness appends observations to `benchmarks/raw_results.jsonl`; make a new copy or rename the old file for an independent benchmark campaign. Do not run benchmarking commands concurrently. A full benchmark temporarily writes an output for every timed run, hashes it after timing, then removes the temporary output. `--keep-method METHOD` retains that method's latest output in the standard results folder. No timing estimates are substituted for unrun full-data cases.

## Folder guide

| Folder/file | Contents |
|---|---|
| `START_HERE.html` | Offline visual overview |
| `REPORT.md` | Measured rankings, ranges, caveats and results |
| `docs/APPROACHES.md` | How each algorithm works, complexity and viva explanation |
| `src/` | Python implementations and plotting code |
| `tests/` | Edge cases, independent reference checks and sampling checks |
| `benchmarks/` | Timing harness, raw logs, summaries, memory and JIT measurements |
| `data/` | Supplied local datasets, metadata and input hashes |
| `results/` | Sorted BED, sampled BEDs, V-plot matrices, expression tables and plots |
| `problem_statement/` | The two supplied PDFs plus extracted text |

The small submission ZIP omits large input archives, raw BED outputs and temporary files. It includes source, documentation, measurements and figures. `.gitignore` excludes large generated data from a future repository. The PDFs mention public GitHub deliverables; this folder is ready to place in a repository, but no public upload has been made.

## Sources

Requirements come from the two PDFs in `problem_statement/`. External technical references, rather than dataset downloads:

* [Python sorting guide](https://docs.python.org/3.13/howto/sorting.html)
* [NumPy sorting behavior](https://numpy.org/doc/stable/reference/generated/numpy.sort.html)
* [SciPy paired Wilcoxon test](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.wilcoxon.html)
* [SciPy independent Mann-Whitney test](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.mannwhitneyu.html)
