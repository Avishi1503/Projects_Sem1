"""Independent edge cases and cross-method checks; no performance assertions."""
import collections
import gzip
import io
import sys
import subprocess
import tempfile
import unittest
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from common import selection, read_bytes, write_order
from kernels import bed_records, bed_offsets_lengths
from sort_bed import sorted_order, external_sort, packed_keys
from rescale import plan_quotas, select_rows, METHODS as SAMPLE_METHODS
from vplot import calculate
from tcga import get_tpm


class Correctness(unittest.TestCase):
    def test_sort_all_methods_custom_order_ties_and_preservation(self):
        lines = [
            b"chr1\t10\t12\tfirst\n",
            b"chr10\t2\t4\textra\tcolumns\n",
            b"chr1\t10\t12\tsecond\n",
            b"chrZ\t0\t9\tremove\n",
            b"chr1\t2\t4\tlast\n",
        ]
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            s = d / "s.tsv"
            s.write_text("chr10\nchr1\n")
            f = d / "a.bed.gz"
            with gzip.open(f, "wb") as stream:
                stream.write(b"".join(lines).rstrip(b"\n"))
            data = read_bytes([f])
            _, names, lens = selection(s)
            rec, seen = bed_records(data, names, lens)
            self.assertEqual(seen, 5)
            expected = b"".join([lines[1], lines[4], lines[0], lines[2]])
            for method in [
                "radix",
                "radix-threaded",
                "radix-jit",
                "radix-jit-threaded",
                "quicksort",
                "heapsort",
                "lexsort",
                "timsort",
                "bucket",
            ]:
                order = sorted_order(rec, method)
                out = d / f"{method}.bed"
                write_order(out, data, rec[:, 3], rec[:, 4], order)
                self.assertEqual(out.read_bytes(), expected, method)
            external_sort([f], s, d / "external.bed", chunk_rows=1)
            self.assertEqual((d / "external.bed").read_bytes(), expected)

    def test_large_coordinates_fallback_and_empty(self):
        rec = np.array(
            [[1, 2**60, 2**60 + 1, 0, 0], [0, 2**60, 2**60 + 2, 0, 0]], np.uint64
        )
        self.assertIsNone(packed_keys(rec))
        self.assertEqual(sorted_order(rec, "radix").tolist(), [1, 0])
        self.assertEqual(len(sorted_order(np.empty((0, 5), np.uint64), "radix")), 0)

    def test_randomized_radix_against_numpy(self):
        from sort_kernels import radix_indices

        rng = np.random.default_rng(213)
        for size in [0, 1, 2, 1000, 10000]:
            for bits in [1, 16, 63]:
                keys = rng.integers(0, 2**bits, size, dtype=np.uint64)
                np.testing.assert_array_equal(radix_indices(keys), np.argsort(keys, kind="stable"))

    def test_parallel_output_matches_serial(self):
        from parallel_writer import write_parallel

        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            f = d / "in.bed"
            f.write_bytes(b"chr1\t9\t10\tlonger payload\nchr1\t1\t2\tx")
            data = read_bytes([f])
            rec = bed_offsets_lengths(data)
            order = np.array([1, 0], np.int64)
            write_order(d / "serial", data, rec[:, 0], rec[:, 1], order)
            unsigned = rec.astype(np.uint64)
            write_parallel(
                d / "parallel", data, unsigned[:, 0], unsigned[:, 1], order, chunk_rows=1
            )
            self.assertEqual((d / "serial").read_bytes(), (d / "parallel").read_bytes())

    def test_vplot_cli_standard_input(self):
        raw = b"2\tchr1\t90\t110\tmotif\t1\t+\tchr1\t98\t103\tread\t5\t+\n"
        with tempfile.TemporaryDirectory() as d:
            out = Path(d) / "matrix.tsv"
            result = subprocess.run(
                [
                    sys.executable,
                    str(Path(__file__).resolve().parents[1] / "src/vplot.py"),
                    "--input",
                    "-",
                    "--output",
                    str(out),
                    "--low",
                    "-2",
                    "--high",
                    "2",
                    "--max-length",
                    "10",
                ],
                input=raw,
                capture_output=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr.decode())
            matrix = np.loadtxt(out, dtype=np.int64)
            self.assertEqual(matrix.shape, (11, 6))
            self.assertEqual(matrix[5, 3], 1)
            self.assertEqual(matrix[:, 1:].sum(), 1)

    def test_invalid_intervals_and_comments(self):
        with self.assertRaises(ValueError):
            bed_offsets_lengths(np.frombuffer(b"chr1\t8\t4\n", np.uint8))
        with self.assertRaises(ValueError):
            bed_offsets_lengths(np.frombuffer(b"chr1\t-1\t4\n", np.uint8))
        rec = bed_offsets_lengths(np.frombuffer(b"# note\ntrack name=x\nchr1\t1\t6\r\n", np.uint8))
        self.assertEqual(rec[:, 2].tolist(), [5])

    def test_sampling_exact_no_replacement_missing_support_and_seed(self):
        with tempfile.TemporaryDirectory() as d:
            ref = Path(d) / "reference.hist"
            ref.write_text("1\t0.2\n2\t0.4\n3\t0.3\n9\t0.1\n")
            lengths = np.array([1] * 20 + [2] * 30 + [3] * 10, dtype=np.int64)
            counts, p, q, info = plan_quotas(np.bincount(lengths), ref, "capped", 0.5)
            self.assertAlmostEqual(info["missing_reference_mass"], 0.1)
            self.assertTrue(np.all(q <= counts))
            for method in SAMPLE_METHODS:
                indices = select_rows(lengths, counts, q, method, 123)
                self.assertEqual(len(set(indices)), len(indices), method)
                np.testing.assert_array_equal(np.bincount(lengths[indices], minlength=len(q)), q)
                np.testing.assert_array_equal(indices, select_rows(lengths, counts, q, method, 123))
            _, _, strict, _ = plan_quotas(np.bincount(lengths), ref, "strict")
            self.assertTrue(np.all(strict <= counts))

    def test_sampling_uniformity_small_population(self):
        # A broad deterministic frequency bound detects first-k selection bias.
        lengths = np.ones(5, dtype=np.int64)
        counts = np.array([0, 5])
        q = np.array([0, 2])
        totals = np.zeros(5)
        for seed in range(1000):
            totals[select_rows(lengths, counts, q, "compiled", seed)] += 1
        self.assertTrue(np.all((totals > 330) & (totals < 470)), totals)

    def test_vplot_known_centers_orientation_and_boundaries(self):
        # region center 100. A half-base fragment center 100.5 bins to 0 on +,
        # and -1 on -. A different chromosome must be excluded.
        rows = [
            b"2\tchr1\t90\t110\tmotif\t1\t+\tchr1\t98\t103\tread\t5\t+",
            b"2\tchr1\t90\t110\tmotif\t1\t-\tchr1\t98\t103\tread\t5\t+",
            b"2\tchr1\t90\t110\tmotif\t1\t+\tchr2\t98\t103\tread\t5\t+",
            b"2\tchr1\t90\t110\tmotif\t1\t+\tchr1\t199\t201\tread\t2\t+",
        ]
        raw = b"\n".join(rows) + b"\n"
        reference = np.zeros((11, 5), np.int64)
        reference[5, 2] = 1
        reference[5, 1] = 1
        for method in ["compiled", "dictionary", "pandas"]:
            matrix, info = calculate(io.BytesIO(raw), method, -2, 2, 10)
            np.testing.assert_array_equal(matrix, reference)
            self.assertEqual(info["input_rows"], 4)
            self.assertEqual(info["counted_rows"], 2)

    def test_tcga_exact_gene_named_tpm_column(self):
        raw = b"# metadata\ngene_id\tgene_name\ttmp\ttpm_unstranded\ng1\tNKX2-10\tx\t999\ng2\tNKX2-1\tx\t53.1718\n"
        for method in ["bytes", "lines", "pandas"]:
            self.assertEqual(get_tpm(raw, "NKX2-1", method), 53.1718)
            with self.assertRaises(ValueError):
                get_tpm(raw, "absent", method)


if __name__ == "__main__":
    unittest.main(verbosity=2)
