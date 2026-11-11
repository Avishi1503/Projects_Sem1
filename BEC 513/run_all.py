"""Run the selected method for each project, then make the plots."""
import argparse
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--projects", type=int, nargs="+", choices=[1, 2, 3, 4], default=[1, 2, 3, 4])
    p.add_argument("--input-dir", type=Path, default=ROOT / "data")
    args = p.parse_args()
    data = args.input_dir.resolve()
    winners_file = ROOT / "benchmarks/winners.json"
    winners = json.loads(winners_file.read_text()) if winners_file.exists() else {}
    defaults = {1: "radix-parallel", 2: "compiled", 3: "compiled", 4: "bytes"}

    def run(script, *extra):
        subprocess.run(
            [sys.executable, str(ROOT / "src" / script), *map(str, extra)], check=True, cwd=ROOT
        )

    for project in args.projects:
        method = winners.get(str(project), {}).get("method", defaults[project])
        print(f"Project {project}: {method}", flush=True)
        if project == 1:
            run(
                "sort_bed.py",
                "--inputs",
                data / "shuf.a.bed.gz",
                data / "shuf.b.bed.gz",
                "--selection",
                data / "selection.tsv",
                "--method",
                method,
            )
        elif project == 2:
            run(
                "rescale.py",
                "--input",
                data / "query.bed.gz",
                "--reference",
                data / "reference.hist",
                "--method",
                method,
            )
            run(
                "rescale.py",
                "--input",
                data / "query.bed.gz",
                "--reference",
                data / "reference.hist",
                "--method",
                method,
                "--policy",
                "strict",
                "--output",
                ROOT / "results/rescale/query.strict.bed",
            )
        elif project == 3:
            run("vplot.py", "--input", data / "mapped.bed.gz", "--method", method)
        else:
            run("tcga.py", "--archive", data / "tcga_source.zip", "--method", method)
        if project in [2, 3, 4]:
            run("make_figures.py", "--project", project)


if __name__ == "__main__":
    main()
