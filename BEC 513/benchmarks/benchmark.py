"""Sequential process benchmarks: wall time, sampled RSS, stage times, output hash.

Do not run competing benchmark processes at the same time. Cold disk cache is
not claimed. Each process starts fresh; pre-existing Numba disk cache is allowed.
"""
import argparse
import csv
import gzip
import hashlib
import json
import os
import platform
import random
import subprocess
import sys
import time
from pathlib import Path
import psutil

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from common import json_write

METHODS = {
    1: [
        "timsort",
        "bucket",
        "quicksort",
        "heapsort",
        "lexsort",
        "radix",
        "radix-threaded",
        "radix-jit",
        "radix-jit-threaded",
        "radix-parallel",
        "external",
    ],
    2: ["python", "reservoir", "priority", "compiled"],
    3: ["dictionary", "pandas", "compiled"],
    4: ["pandas", "lines", "bytes"],
}


def prepare_samples():
    directory = ROOT / "data/screening"
    directory.mkdir(exist_ok=True)
    for name in ["shuf.a.bed.gz", "shuf.b.bed.gz", "query.bed.gz"]:
        target = directory / name
        if not target.exists():
            with gzip.open(ROOT / "data" / name, "rb") as f, gzip.open(
                target, "wb", compresslevel=1
            ) as g:
                for _ in range(100000 if name.startswith("shuf") else 200000):
                    line = f.readline()
                    if not line:
                        break
                    g.write(line)
    return directory


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def run_job(project, method, scope, repeat, keep=False):
    temp = ROOT / "tmp/benchmark"
    temp.mkdir(parents=True, exist_ok=True)
    base = f"p{project}.{scope}.{method}.{repeat}"
    output = temp / (base + (".bed" if project < 3 else ".tsv"))
    metric = temp / (base + ".metrics.json")
    log = ROOT / "benchmarks/logs" / f"{base}.log"
    log.parent.mkdir(exist_ok=True)
    script = {1: "sort_bed.py", 2: "rescale.py", 3: "vplot.py", 4: "tcga.py"}[project]
    cmd = [
        sys.executable,
        str(ROOT / "src" / script),
        "--method",
        method,
        "--output",
        str(output),
        "--metrics",
        str(metric),
    ]
    if project == 1 and scope == "screening":
        data = prepare_samples()
        cmd += ["--inputs", str(data / "shuf.a.bed.gz"), str(data / "shuf.b.bed.gz")]
    if project == 2 and scope == "screening":
        cmd += ["--input", str(prepare_samples() / "query.bed.gz")]
    if project == 3:
        cmd += ["--input", str(ROOT / "data/mapped.bed.gz")]
        if scope == "screening":
            cmd += ["--limit", "200000"]
    if project == 4 and scope == "screening":
        cmd += ["--limit", "12"]
    env = dict(os.environ)
    env.update(
        OMP_NUM_THREADS="1",
        OPENBLAS_NUM_THREADS="1",
        MKL_NUM_THREADS="1",
        NUMEXPR_NUM_THREADS="1",
        NUMBA_NUM_THREADS="1",
        PYTHONHASHSEED="0",
    )
    if project == 1 and method == "radix-parallel":
        env["NUMBA_NUM_THREADS"] = "4"
    before = psutil.virtual_memory().available
    started = time.perf_counter()
    peak = 0
    peak_private = 0
    with log.open("w", encoding="utf-8") as f:
        proc = subprocess.Popen(cmd, stdout=f, stderr=subprocess.STDOUT, env=env)
        process = psutil.Process(proc.pid)
        while proc.poll() is None:
            try:
                mem = process.memory_info()
                peak = max(peak, mem.rss)
                peak_private = max(peak_private, getattr(mem, "private", 0))
            except psutil.NoSuchProcess:
                pass
            time.sleep(0.05)
    wall = time.perf_counter() - started
    row = dict(
        project=project,
        method=method,
        scope=scope,
        repeat=repeat,
        wall_seconds=wall,
        peak_rss_mb=peak / 2**20,
        peak_private_mb=peak_private / 2**20,
        available_before_mb=before / 2**20,
        exit_code=proc.returncode,
        command=cmd,
        log=str(log.relative_to(ROOT)),
    )
    if proc.returncode:
        row["error"] = log.read_text()[-6000:]
    else:
        row.update(json.loads(metric.read_text()))
        row["output_bytes"] = output.stat().st_size
        row["output_sha256"] = sha256(output)
        if project == 3:
            row["matrix_sha256"] = sha256(output.with_suffix(".npy"))
        if keep:
            final = {
                1: ROOT / "results/sorted_bed_file_per_sample/X.provided.sorted.bed",
                2: ROOT / "results/rescale/query.rescaled.bed",
                3: ROOT / "results/vplot/vplot.matrix.tsv",
                4: ROOT / "results/tcga/NKX2-1.tsv",
            }[project]
            final.parent.mkdir(parents=True, exist_ok=True)
            os.replace(output, final)
            for suffix in [".npy", ".summary.json", ".hist.tsv"]:
                other = output.with_suffix(suffix)
                if other.exists():
                    os.replace(other, final.with_suffix(suffix))
            row["retained_output"] = str(final.relative_to(ROOT))
    for p in temp.glob(base + "*"):
        if p.is_file():
            p.unlink()
    return row


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--project", type=int, choices=[1, 2, 3, 4], required=True)
    p.add_argument("--scope", choices=["screening", "full"], default="screening")
    p.add_argument("--methods", nargs="+")
    p.add_argument("--repeats", type=int, default=3)
    p.add_argument("--keep-method")
    args = p.parse_args()
    methods = args.methods or METHODS[args.project]
    if any(m not in METHODS[args.project] for m in methods):
        p.error("Invalid method")
    result_file = ROOT / "benchmarks/raw_results.jsonl"
    jobs = [(m, r) for r in range(1, args.repeats + 1) for m in methods]
    random.Random(20261105 + args.project).shuffle(jobs)
    for m, r in jobs:
        print(f"START project={args.project} scope={args.scope} method={m} repeat={r}", flush=True)
        row = run_job(args.project, m, args.scope, r, m == args.keep_method)
        with result_file.open("a", encoding="utf-8") as f:
            f.write(json.dumps(row, allow_nan=False) + "\n")
        print(
            json.dumps(
                {
                    k: row.get(k)
                    for k in [
                        "project",
                        "method",
                        "repeat",
                        "wall_seconds",
                        "peak_rss_mb",
                        "exit_code",
                        "error",
                    ]
                }
            ),
            flush=True,
        )
        if row["exit_code"]:
            raise SystemExit("Benchmark failed; inspect its log")


if __name__ == "__main__":
    main()
