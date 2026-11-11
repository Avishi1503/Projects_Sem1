"""Sample BED rows to match a reference fragment-length distribution."""
import argparse
import random
import time
from pathlib import Path
import numpy as np
from common import ROOT, read_bytes, write_order, emit_metrics, json_write
from kernels import bed_offsets_lengths, exact_select

METHODS = ["compiled", "python", "priority", "reservoir"]


def plan_quotas(counts, reference, policy="capped", fraction=0.30):
    """Calculate how many rows to keep at each fragment length.

    Strict scales all shared lengths together. Capped uses fraction*N as
    the target, then limits each quota to the rows available.
    """
    refs = np.atleast_2d(np.loadtxt(reference, dtype=float))
    if refs.shape[1] != 2 or not np.all(np.isfinite(refs)):
        raise ValueError("Reference must contain two finite numeric columns")
    if (
        np.any(refs[:, 0] < 0)
        or np.any(refs[:, 0] != np.floor(refs[:, 0]))
        or np.any(refs[:, 1] < 0)
    ):
        raise ValueError("Reference lengths must be integers and weights nonnegative")
    size = max(len(counts), int(refs[:, 0].max()) + 1)
    if size > 10_000_000:
        raise ValueError("Length range too large for dense histogram")
    counts = np.pad(counts, (0, size - len(counts)))
    p = np.zeros(size)
    np.add.at(p, refs[:, 0].astype(int), refs[:, 1])
    if p.sum() <= 0:
        raise ValueError("Reference mass is zero")
    p /= p.sum()
    support = (counts > 0) & (p > 0)
    missing = float(p[(counts == 0) & (p > 0)].sum())
    if not np.any(support):
        raise ValueError("Query and reference have no common length support")
    if policy == "strict":
        # p is left on the original scale: alpha*p gives the same quotas as
        # renormalizing on common support and rescaling alpha appropriately.
        scale = float(np.min(counts[support] / p[support]))
    else:
        scale = float(counts.sum() * fraction)
    quotas = np.minimum(counts, np.floor(scale * p + 1e-10).astype(np.int64))
    if quotas.sum() == 0:
        raise ValueError("No rows survive integer quotas; use capped policy or a larger dataset")
    realized = quotas / quotas.sum()
    info = dict(
        policy=policy,
        requested_scale=scale,
        query_rows=int(counts.sum()),
        selected_rows=int(quotas.sum()),
        missing_reference_mass=missing,
        total_variation=float(0.5 * np.abs(realized - p).sum()),
        saturated_bins=int(np.sum((scale * p > counts) & (p > 0))),
        integer_rounding="floor",
        replacement=False,
    )
    return counts, p, quotas, info


def select_rows(lengths, counts, quotas, method, seed):
    if method == "compiled":
        return exact_select(lengths, counts, quotas, seed)
    rng = np.random.default_rng(seed)
    if method == "priority":
        # Independent random priorities give uniform subsets within each bin.
        order = np.lexsort((rng.random(len(lengths)), lengths))
        cumulative = np.concatenate(([0], np.cumsum(counts)))
        chosen = [
            order[int(cumulative[k]) : int(cumulative[k] + q)] for k, q in enumerate(quotas) if q
        ]
        result = np.concatenate(chosen) if chosen else np.empty(0, np.int64)
        result.sort()
        return result
    if method == "reservoir":
        buckets = [np.empty(int(q), np.int64) for q in quotas]
        seen = np.zeros(len(counts), np.int64)
        for i, k in enumerate(lengths):
            k = int(k)
            q = int(quotas[k])
            n = int(seen[k])
            seen[k] += 1
            if n < q:
                buckets[k][n] = i
            elif q:
                j = int(rng.integers(n + 1))
                if j < q:
                    buckets[k][j] = i
        result = np.concatenate(buckets)
        result.sort()
        return result
    generator = random.Random(seed)
    remaining = counts.copy()
    need = quotas.copy()
    chosen = []
    for i, k in enumerate(lengths):
        k = int(k)
        if need[k] > 0 and (need[k] == remaining[k] or generator.random() * remaining[k] < need[k]):
            chosen.append(i)
            need[k] -= 1
        remaining[k] -= 1
    return np.array(chosen, np.int64)


def run(args):
    started = time.perf_counter()
    t = time.perf_counter()
    data = read_bytes([args.input], args.limit)
    read_s = time.perf_counter() - t
    t = time.perf_counter()
    rec = bed_offsets_lengths(data)
    parse_s = time.perf_counter() - t
    lengths = np.ascontiguousarray(rec[:, 2])
    if len(lengths) == 0 or lengths.max() > 10_000_000:
        raise ValueError("Empty input or fragment length too large for dense histogram")
    t = time.perf_counter()
    counts = np.bincount(lengths)
    counts, p, quotas, info = plan_quotas(counts, args.reference, args.policy, args.fraction)
    plan_s = time.perf_counter() - t
    t = time.perf_counter()
    chosen = select_rows(lengths, counts, quotas, args.method, args.seed)
    sample_s = time.perf_counter() - t
    # Check exact quotas and distinct original records outside the selection timer.
    actual = np.bincount(lengths[chosen], minlength=len(counts))
    if not np.array_equal(actual, quotas) or (
        len(chosen) > 1 and np.any(chosen[1:] <= chosen[:-1])
    ):
        raise AssertionError("Invalid sample")
    t = time.perf_counter()
    write_order(args.output, data, rec[:, 0], rec[:, 1], chosen)
    write_s = time.perf_counter() - t
    hist_path = Path(args.output).with_suffix(".hist.tsv")
    hist = np.column_stack((np.arange(len(counts)), counts, p, quotas, actual / actual.sum()))
    np.savetxt(
        hist_path,
        hist,
        delimiter="\t",
        fmt=["%d", "%d", "%.12g", "%d", "%.12g"],
        header="length\tquery_count\treference_probability\tsampled_count\tsampled_probability",
        comments="",
    )
    info.update(
        project=2,
        method=args.method,
        seed=args.seed,
        requested_fraction=args.fraction if args.policy == "capped" else None,
        read_seconds=read_s,
        parse_seconds=parse_s,
        plan_seconds=plan_s,
        sample_seconds=sample_s,
        write_seconds=write_s,
        quota_check=True,
        distinct_indices=True,
    )
    json_write(Path(args.output).with_suffix(".summary.json"), info)
    emit_metrics(args, info, started)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", default=str(ROOT / "data/query.bed.gz"))
    p.add_argument("--reference", default=str(ROOT / "data/reference.hist"))
    p.add_argument("--output", default=str(ROOT / "results/rescale/query.rescaled.bed"))
    p.add_argument("--method", choices=METHODS, default="compiled")
    p.add_argument("--policy", choices=["strict", "capped"], default="capped")
    p.add_argument("--fraction", type=float, default=0.30)
    p.add_argument("--seed", type=int, default=2026)
    p.add_argument("--limit", type=int, default=0)
    p.add_argument("--metrics")
    args = p.parse_args()
    if not 0 < args.fraction <= 1 or args.limit < 0 or not 0 <= args.seed < 2**32:
        p.error("Require 0 < fraction <= 1, limit >= 0, and uint32 seed")
    run(args)
