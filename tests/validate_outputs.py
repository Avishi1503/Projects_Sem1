"""Independent checks on delivered files, outside benchmark timing."""
import collections
import csv
import hashlib
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main():
    checks = {}
    observations = [
        json.loads(line)
        for line in (ROOT / "benchmarks/raw_results.jsonl").read_text().splitlines()
        if line.strip()
    ]
    hashes = {
        r["output_sha256"]
        for r in observations
        if r["project"] == 1 and r["scope"] == "full" and r["exit_code"] == 0
    }
    actual = sha(ROOT / "results/sorted_bed_file_per_sample/X.provided.sorted.bed")
    assert hashes == {actual}, "Sorted outputs disagree"
    checks["sort"] = dict(
        sha256=actual,
        all_full_methods_equal=True,
        independent_external_merge_included=any(
            r["project"] == 1 and r["method"] == "external" and r["scope"] == "full"
            for r in observations
        ),
    )
    checks["rescale"] = {}
    for stem in ["query.rescaled", "query.strict"]:
        file = ROOT / "results/rescale" / f"{stem}.bed"
        histogram = collections.Counter()
        with file.open("rb") as f:
            for line in f:
                fields = line.split(b"\t", 3)
                histogram[int(fields[2]) - int(fields[1])] += 1
        with file.with_suffix(".hist.tsv").open() as f:
            expected = {
                int(r["length"]): int(r["sampled_count"]) for r in csv.DictReader(f, delimiter="\t")
            }
        assert dict(histogram) == {k: v for k, v in expected.items() if v}, stem
        checks["rescale"][stem] = dict(
            output_rows=sum(histogram.values()), written_histogram_verified=True, sha256=sha(file)
        )
    matrix_file = ROOT / "results/vplot/vplot.matrix.tsv"
    wide = np.loadtxt(matrix_file, dtype=np.int64)
    binary = np.load(matrix_file.with_suffix(".npy"))
    np.testing.assert_array_equal(wide[:, 0], np.arange(len(binary)))
    np.testing.assert_array_equal(wide[:, 1:], binary)
    hashes = {
        r["output_sha256"]
        for r in observations
        if r["project"] == 3 and r["scope"] == "full" and r["exit_code"] == 0
    }
    assert hashes == {sha(matrix_file)}
    checks["vplot"] = dict(
        shape=list(binary.shape),
        matrix_sum=int(binary.sum()),
        wide_equals_binary=True,
        all_full_methods_equal=True,
    )
    table = ROOT / "results/tcga/NKX2-1.tsv"
    hashes = {
        r["output_sha256"]
        for r in observations
        if r["project"] == 4 and r["scope"] == "full" and r["exit_code"] == 0
    }
    assert hashes == {sha(table)}
    with table.open() as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    hint = [r for r in rows if r["file_id"] == "e6e83fd4-2f80-44e1-8dce-588a37308dbb"]
    assert len(hint) == 1 and float(hint[0]["tpm"]) == 53.1718
    assert len({r["file_id"] for r in rows}) == len(rows)
    checks["tcga"] = dict(
        rows=len(rows), all_full_methods_equal=True, unique_file_ids=True, hint_tpm=53.1718
    )
    (ROOT / "results/validation.json").write_text(json.dumps(checks, indent=2))
    print(json.dumps(checks, indent=2))


if __name__ == "__main__":
    main()
