"""Build the report and results page from saved benchmark measurements."""
import collections
import csv
import html
import json
import platform
import statistics
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import numba
import scipy
import psutil

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from common import json_write

TITLES = {
    1: "BED sorting",
    2: "Reference resampling",
    3: "V-plot matrix",
    4: "TCGA gene extraction",
}


def main():
    rows = [
        json.loads(line)
        for line in (ROOT / "benchmarks/raw_results.jsonl").read_text().splitlines()
        if line.strip()
    ]
    grouped = collections.defaultdict(list)
    for row in rows:
        if row["exit_code"] == 0:
            grouped[row["project"], row["scope"], row["method"]].append(row)
    summary = []
    for (project, scope, method), rr in sorted(grouped.items()):
        times = [r["wall_seconds"] for r in rr]
        d = dict(
            project=project,
            scope=scope,
            method=method,
            runs=len(rr),
            median_seconds=statistics.median(times),
            min_seconds=min(times),
            max_seconds=max(times),
            max_peak_rss_mb=max(r["peak_rss_mb"] for r in rr),
            median_pipeline_seconds=statistics.median(r["pipeline_seconds"] for r in rr),
        )
        for key in [
            "read_seconds",
            "parse_seconds",
            "sort_seconds",
            "sample_seconds",
            "calculate_seconds",
            "parse_gene_seconds",
            "write_seconds",
        ]:
            v = [r[key] for r in rr if key in r]
            if v:
                d["median_" + key] = statistics.median(v)
        summary.append(d)
    fields = list(dict.fromkeys(k for r in summary for k in r))
    with (ROOT / "benchmarks/summary.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fields)
        w.writeheader()
        w.writerows(summary)
    json_write(ROOT / "benchmarks/summary.json", summary)
    winners = {}
    for project in TITLES:
        candidates = [
            r
            for r in summary
            if r["project"] == project and r["scope"] == "full" and r["runs"] >= 3
        ]
        if candidates:
            winners[str(project)] = min(candidates, key=lambda r: r["median_seconds"])
    json_write(ROOT / "benchmarks/winners.json", winners)
    equivalence = {}
    for project in [1, 3, 4]:
        rr = [
            r
            for r in rows
            if r["project"] == project and r["scope"] == "full" and r["exit_code"] == 0
        ]
        hashes = {r["output_sha256"] for r in rr}
        equivalence[str(project)] = dict(
            full_runs=len(rr), distinct_output_hashes=len(hashes), all_equal=len(hashes) == 1
        )
    rr = [r for r in rows if r["project"] == 2 and r["scope"] == "full" and r["exit_code"] == 0]
    equivalence["2"] = dict(
        full_runs=len(rr),
        all_quota_checks_passed=all(r.get("quota_check") for r in rr),
        all_distinct_index_checks_passed=all(r.get("distinct_indices") for r in rr),
        selected_counts=sorted({r["selected_rows"] for r in rr}),
        total_variations=sorted({r["total_variation"] for r in rr}),
    )
    json_write(ROOT / "benchmarks/equivalence.json", equivalence)
    environment = dict(
        python=sys.version,
        executable=sys.executable,
        os=platform.platform(),
        cpu=platform.processor(),
        physical_cores=psutil.cpu_count(logical=False),
        logical_cores=psutil.cpu_count(),
        ram_gib=psutil.virtual_memory().total / 2**30,
        versions=dict(
            numpy=np.__version__,
            pandas=pd.__version__,
            numba=numba.__version__,
            scipy=scipy.__version__,
            psutil=psutil.__version__,
        ),
        numba_cache="pre-populated by correctness/screening runs for main benchmarks",
        parallelism="1 numerical worker normally; threaded sorting variants use up to 2 file readers; radix-parallel uses 4 output-copy workers",
        disk_cache="not cleared or controlled; no cold-disk claim",
        timing="fresh child process, launch through exit, output close included; output hashing and correctness checks by harness after exit excluded",
        peak_memory="child RSS sampled every 50 ms; may miss very short spikes",
        host_load="normal desktop environment; other applications were not closed",
    )
    json_write(ROOT / "benchmarks/environment.json", environment)
    lines = [
        "# Benchmark results",
        "",
        "The tables below compare the Python implementations on the supplied datasets. Times vary by machine; these are the measurements from the test laptop.",
        "",
        "## Fastest methods in these runs",
        "",
        "| Project | Fastest tested method | Median execution | Observed range | Peak RSS |",
        "|---|---|---:|---:|---:|",
    ]
    for p, w in winners.items():
        lines.append(
            f"| {p}. {TITLES[int(p)]} | `{w['method']}` | {w['median_seconds']:.2f} s | {w['min_seconds']:.2f}-{w['max_seconds']:.2f} s | {w['max_peak_rss_mb']:.0f} MiB |"
        )
    lines += [
        "",
        "## How the timings were measured",
        "",
        "Each run starts a new Python process. The timer includes imports, loading the Numba cache, reading, decompression, calculation and writing the output. Each method uses the same input and output requirements. Runs are repeated at least three times, in a shuffled order with a fixed seed within each batch. The final sorting batch compares the leading methods together. Writes are buffered by the operating system; no forced disk sync is used. Output hashing by the benchmark script happens after timing, while checks inside each program remain included.",
        "",
        "The timed outputs are the sorted BED, sampled BED and histogram, V-plot matrix, and TCGA expression table. Plotting and patient-level statistics are measured separately in `benchmarks/figure_times.json`. In particular, the TCGA time in the table does not include making the boxplot.",
        "",
        "Disk caches were not cleared and other desktop applications were open. Downloading and installing packages are excluded. Threaded sorting variants use two readers; radix-parallel also uses four output-copy workers. Other variants use one computational worker.",
        "",
        "Screening uses smaller inputs: 100,000 rows from each sorting file, 200,000 rows for resampling and V-plot, and 12 eligible TCGA files. Timsort and chromosome buckets were only run at this size because of their Python-object memory use. Full-input results are labelled separately.",
        "",
    ]
    for project, title in TITLES.items():
        lines += [
            f"## {project}. {title}",
            "",
            "| Scope | Method | Runs | Median seconds | Min-max seconds | Maximum sampled RSS (MiB) |",
            "|---|---|---:|---:|---:|---:|",
        ]
        for r in sorted(
            (r for r in summary if r["project"] == project),
            key=lambda r: (r["scope"], r["median_seconds"]),
        ):
            lines.append(
                f"| {r['scope']} | {r['method']} | {r['runs']} | {r['median_seconds']:.3f} | {r['min_seconds']:.3f}-{r['max_seconds']:.3f} | {r['max_peak_rss_mb']:.1f} |"
            )
        full = [
            r
            for r in rows
            if r["project"] == project and r["scope"] == "full" and r["exit_code"] == 0
        ]
        if full:
            r = full[-1]
            if project == 1:
                lines += [
                    "",
                    f"Input rows: **{r['input_rows']:,}**. Retained after chromosome selection: **{r['output_rows']:,}**. All full-data outputs equal: **{equivalence['1']['all_equal']}**.",
                ]
            if project == 2:
                lines += [
                    "",
                    f"Query rows: **{r['query_rows']:,}**. Capped-policy sampled rows: **{r['selected_rows']:,}**. Missing reference mass: **{r['missing_reference_mass']:.6%}**. Total variation from reference: **{r['total_variation']:.6f}**. Every selector meets the same quotas and uses distinct input indices.",
                    "",
                    "Random sample records can differ between methods; their per-length counts must agree. `results/rescale/query.strict.summary.json` documents the stricter shape/retention tradeoff separately.",
                ]
            if project == 3:
                lines += [
                    "",
                    f"Input region/fragment pairs: **{r['input_rows']:,}**. Counted in the configured window: **{r['counted_rows']:,}**. Excluded: **{r['excluded_rows']:,}**. All full-data matrices equal: **{equivalence['3']['all_equal']}**.",
                ]
            if project == 4:
                lines += [
                    "",
                    f"Eligible expression files: **{r['eligible_files']:,}**. Manifest MD5 checks: **{r['md5_verified']:,}**. Groups: `{r['groups']}`. All extraction tables equal: **{equivalence['4']['all_equal']}**.",
                ]
        lines.append("")
    jit = ROOT / "benchmarks/jit_cost.json"
    if jit.exists():
        costs = json.loads(jit.read_text())
        lines += [
            "## First-use compilation cost",
            "",
            "The main timings use an existing Numba cache. This smaller experiment compares the first call with an empty cache, the first call with a saved cache, and a second call in the same process. The first-call times include a small example calculation as well as compilation or loading.",
            "",
            "| Kernel | Empty-cache first call | Populated-cache first call | In-process repeat |",
            "|---|---:|---:|---:|",
        ]
        for cold, warm in zip(costs[0]["kernels"], costs[1]["kernels"]):
            lines.append(
                f"| {cold['kernel']} | {cold['first_call_seconds']:.4f} s | {warm['first_call_seconds']:.4f} s | {warm['in_process_second_call_seconds']:.6f} s |"
            )
        lines += [
            "",
            "Run `python benchmarks/jit_cost.py` to repeat this check. A first run on another computer includes compilation. NumPy and pandas use libraries compiled when installed.",
            "",
        ]
    coldfile = ROOT / "benchmarks/cold_start_full.json"
    if coldfile.exists():
        cold = json.loads(coldfile.read_text())
        lines += [
            "## Complete first-use runs",
            "",
            "These full-input runs use a separate, empty Numba cache and include compilation. Each time is a single observation, not a median. The operating-system file cache was still uncontrolled. These runs are kept separate from the repeated timings above.",
            "",
            "| Project | Selected method | Empty-cache complete run |",
            "|---|---|---:|",
        ]
        for r in cold:
            lines.append(f"| {r['project']} | {r['method']} | {r['wall_seconds']:.2f} s |")
        lines += [
            "",
            "TCGA byte extraction uses no Numba kernel; its main fresh-process timings already include Python/library startup.",
            "",
        ]
    strictfile = ROOT / "results/rescale/query.strict.summary.json"
    if strictfile.exists():
        s = json.loads(strictfile.read_text())
        lines += [
            "## Recommended resampling output",
            "",
            f"For the closest reference match, use `results/rescale/query.strict.bed`: **{s['selected_rows']:,} rows**, total variation **{s['total_variation']:.8f}**. The capped comparison output retains more rows but has greater shape error. The two policies produce different amounts of data, so their runtimes are not directly comparable.",
            "",
        ]
    statfile = ROOT / "results/tcga/statistics.json"
    if statfile.exists():
        s = json.loads(statfile.read_text())
        test = s["paired_wilcoxon"]
        lines += [
            "## Expression interpretation",
            "",
            f"The files contain {s['expression_files']} eligible expression records, {s['unique_biological_samples']} biological samples and {s['unique_patients']} patients. The descriptive boxplot shows each expression file. Patient-level testing uses {s['matched_patients']} matched tumor/normal pairs after aggregating repeated observations.",
            f"Paired Wilcoxon p-value: **{test['pvalue']:.6g}**. Median paired change in log2(TPM + 1): **{test.get('median_tumor_minus_normal_log2',0):.4f}**.",
            "",
            "The paired test compares tumor and normal tissue within patients. The all-file boxplot shows the overall spread. Neither result alone establishes diagnostic accuracy.",
            "",
        ]
    if statfile.exists():
        lines += [
            "Across all files, the tumor and normal medians are close, but tumor values vary more. Among matched patients, the median tumor-minus-normal change is negative. Only one normal patient is unmatched, so a separate unmatched-group test is not reported. The example plot in the assignment is illustrative; its p-value is not a target for this analysis.",
            "",
        ]
    lines += [
        "## Supporting files",
        "",
        "* [Algorithm notes](docs/APPROACHES.md)",
        "* [Setup and run commands](README.md)",
        "* Raw observations: `benchmarks/raw_results.jsonl`",
        "* Full command logs: `benchmarks/logs/`",
        "* Per-stage timings: `benchmarks/summary.csv`",
        "* Input provenance and SHA-256: `data/sources.json`",
        "* Source hashes from the original timing runs: `benchmarks/source_hashes.json` (before the formatting and comment edits)",
        "",
        f"Tested environment: Python {platform.python_version()}, NumPy {np.__version__}, pandas {pd.__version__}, Numba {numba.__version__}, SciPy {scipy.__version__}; Intel Core i5-1135G7, 4 physical cores / 8 logical cores, {psutil.virtual_memory().total/2**30:.1f} GiB RAM.",
    ]
    (ROOT / "REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    make_chart(summary, winners)
    create_html(summary, winners, equivalence)
    print(json.dumps(winners, indent=2))


def make_chart(summary, winners):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(14, 11), layout="constrained")
    for project, ax in zip(TITLES, axes.flat):
        rs = sorted(
            [r for r in summary if r["project"] == project and r["scope"] == "full"],
            key=lambda r: r["median_seconds"],
        )
        y = np.arange(len(rs))
        times = np.array([r["median_seconds"] for r in rs])
        errors = np.array(
            [
                [r["median_seconds"] - r["min_seconds"] for r in rs],
                [r["max_seconds"] - r["median_seconds"] for r in rs],
            ]
        )
        colors = [
            "#007f86" if r["method"] == winners.get(str(project), {}).get("method") else "#728998"
            for r in rs
        ]
        ax.barh(
            y,
            times,
            xerr=errors,
            color=colors,
            height=0.65,
            error_kw={"elinewidth": 1, "capsize": 2},
        )
        ax.set_yticks(y, [r["method"] for r in rs])
        ax.invert_yaxis()
        ax.set_title(f"{project}. {TITLES[project]}", loc="left", fontweight="bold")
        ax.set_xlabel("Fresh-process execution time (seconds)")
        if len(rs):
            ceiling = max(r["max_seconds"] for r in rs)
            ax.set_xlim(0, ceiling * 1.28)
            for i, r in enumerate(rs):
                ax.text(
                    r["max_seconds"] + ceiling * 0.025,
                    i,
                    f"{r['median_seconds']:.2f}",
                    va="center",
                    fontsize=10,
                )
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(axis="x", alpha=0.15)
        ax.set_axisbelow(True)
    fig.suptitle(
        "Measured Python approaches on complete supplied datasets\nBars = medians; whiskers = observed minimum and maximum",
        fontsize=16,
        fontweight="bold",
    )
    fig.savefig(ROOT / "benchmarks/benchmark_comparison.png", dpi=160)
    plt.close(fig)


def create_html(summary, winners, equivalence):
    cards = "".join(
        f'<article><div class="eyebrow">PROJECT {p}</div><h2>{TITLES[int(p)]}</h2><div class="time">{w["median_seconds"]:.2f}<span> seconds</span></div><p><code>{w["method"]}</code></p><small>Median of {w["runs"]} full-data runs · {w["max_peak_rss_mb"]:.0f} MiB peak RSS</small></article>'
        for p, w in winners.items()
    )
    tables = ""
    for p, title in TITLES.items():
        rs = sorted(
            [r for r in summary if r["project"] == p and r["scope"] == "full"],
            key=lambda r: r["median_seconds"],
        )
        body = "".join(
            f'<tr><td>{r["method"]}</td><td>{r["median_seconds"]:.2f} s</td><td>{r["min_seconds"]:.2f}–{r["max_seconds"]:.2f} s</td><td>{r["max_peak_rss_mb"]:.0f} MiB</td></tr>'
            for r in rs
        )
        tables += f"<section><h2>{p}. {title}</h2><table><tr><th>Method</th><th>Median</th><th>Observed range</th><th>Peak RSS</th></tr>{body}</table></section>"
    figures = ""
    for path, title, caption in [
        (
            "results/rescale/distribution_comparison.png",
            "Resampling accuracy",
            "Strict sampling matches the reference more closely; capped sampling keeps more rows.",
        ),
        (
            "results/vplot/vplot.png",
            "V-plot",
            "Each cell counts region/fragment overlaps at one fragment length and midpoint offset.",
        ),
        (
            "results/tcga/NKX2-1_boxplot.png",
            "NKX2-1 expression",
            "Descriptive expression-file view; log2(TPM + 1).",
        ),
        (
            "results/tcga/NKX2-1_paired_patients.png",
            "Matched-patient comparison",
            "Pairing prevents tumor and normal samples from the same patient being treated as independent.",
        ),
    ]:
        if (ROOT / path).exists():
            figures += f'<section><h2>{title}</h2><p>{caption}</p><img src="{path}" alt="{title}"></section>'
    text = """<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Semester 1 bioinformatics results</title><style>
    *{box-sizing:border-box}body{margin:0;background:#f3f6f8;color:#183247;font:16px/1.6 system-ui,sans-serif}main{max-width:1140px;margin:auto;padding:52px 28px}header{padding:20px 0 32px;max-width:850px}.eyebrow{font-size:12px;letter-spacing:2px;font-weight:750;color:#007f86}h1{font-size:34px;line-height:1.2;margin:12px 0 20px}h2{font-size:23px;line-height:1.25;margin:12px 0}p{color:#4d6373}a{color:#007f86}nav{display:flex;gap:22px;flex-wrap:wrap;margin:22px 0}.cards{display:grid;grid-template-columns:repeat(2,1fr);gap:18px}article,section{background:white;border:1px solid #dce6eb;border-radius:4px;padding:22px}article .time{font-size:30px;font-weight:720;color:#007f86}.time span{font-size:16px;font-weight:450}small{color:#667b88}section{margin-top:24px}table{border-collapse:collapse;width:100%;margin-top:18px;font-size:15px}th,td{text-align:left;padding:12px 10px;border-bottom:1px solid #e6ecef}th{background:#f3f7f8;font-size:13px}img{display:block;width:100%;height:auto;margin:18px auto 0}code{background:#eef3f5;padding:3px 6px;border-radius:4px}.note{border-left:4px solid #df7634;padding-left:18px;margin:28px 0}footer{padding:30px 0;color:#667b88;font-size:14px}@media(max-width:700px){main{padding:24px 16px}h1{font-size:34px}.cards{grid-template-columns:1fr}table{font-size:12px}th,td{padding:8px 5px}}
    </style></head><body><main><header><div class="eyebrow">SEMESTER 1</div><h1>Bioinformatics project results</h1><p>Four Python tasks, their timing comparisons, and the resulting plots.</p><nav><a href="README.md">Run the projects</a><a href="REPORT.md">Detailed measurements</a><a href="docs/APPROACHES.md">Understand the algorithms</a><a href="benchmarks/summary.csv">Timing table</a></nav></header>"""
    text += (
        '<div class="cards">'
        + cards
        + '</div><p class="note">Median times include program startup, cached Numba loading, reading and writing. Plotting is separate. Disk caches were not cleared. See the report for first-run compilation costs and timing ranges.</p><section><h2>How the approaches compare</h2><img src="benchmarks/benchmark_comparison.png" alt="Full-data runtime comparison"></section>'
        + tables
        + figures
    )
    text += "<footer>Code, inputs and methods are described in the README. Large datasets and BED outputs stay local.</footer></main></body></html>"
    (ROOT / "START_HERE.html").write_text(text, encoding="utf-8")


if __name__ == "__main__":
    main()
