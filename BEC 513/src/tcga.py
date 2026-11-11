"""Extract a gene's TPM from tumor and normal samples in a TCGA archive."""
import argparse
import csv
import hashlib
import io
import math
import tarfile
import time
import zipfile
from pathlib import Path, PurePosixPath
from common import ROOT, json_write, emit_metrics

METHODS = ["bytes", "lines", "pandas"]
GROUPS = {"Primary Tumor", "Solid Tissue Normal"}


def get_tpm(raw, gene, method):
    token = gene.encode("utf-8")
    header = (
        next((line for line in io.BytesIO(raw) if not line.startswith(b"#")), b"")
        .rstrip(b"\r\n")
        .split(b"\t")
    )
    name_col = header.index(b"gene_name")
    tpm_col = header.index(b"tpm_unstranded")
    if method == "pandas":
        import pandas as pd

        frame = pd.read_csv(
            io.BytesIO(raw),
            sep="\t",
            comment="#",
            usecols=["gene_name", "tpm_unstranded"],
            float_precision="round_trip",
        )
        hits = frame.loc[frame["gene_name"] == gene, "tpm_unstranded"]
        if len(hits) != 1:
            raise ValueError(f"{gene}: expected one row, found {len(hits)}")
        value = float(hits.iloc[0])
    elif method == "lines":
        values = []
        for line in io.BytesIO(raw):
            fields = line.rstrip(b"\r\n").split(b"\t")
            if len(fields) > max(name_col, tpm_col) and fields[name_col] == token:
                values.append(float(fields[tpm_col]))
        if len(values) != 1:
            raise ValueError(f"{gene}: expected one row, found {len(values)}")
        value = values[0]
    else:
        needle = b"\t" + token + b"\t"
        at = 0
        values = []
        while True:
            match = raw.find(needle, at)
            if match < 0:
                break
            start = raw.rfind(b"\n", 0, match) + 1
            end = raw.find(b"\n", match)
            if end < 0:
                end = len(raw)
            fields = raw[start:end].rstrip(b"\r").split(b"\t")
            if len(fields) > max(name_col, tpm_col) and fields[name_col] == token:
                values.append(float(fields[tpm_col]))
            at = end + 1
        if len(values) != 1:
            raise ValueError(f"{gene}: expected one row, found {len(values)}")
        value = values[0]
    if not math.isfinite(value) or value < 0:
        raise ValueError("TPM must be finite and nonnegative")
    return value


def extract(archive, gene, method, limit=0, verify=True):
    rows = []
    parse_s = 0.0
    hash_s = 0.0
    counts = {
        "expression_files_in_archive": 0,
        "eligible_files": 0,
        "excluded_expression_files": 0,
        "md5_verified": 0,
    }
    with zipfile.ZipFile(archive) as z:
        with z.open("gdc_sample_sheet.tsv") as stream:
            metadata = {}
            for row in csv.DictReader(
                io.TextIOWrapper(stream, encoding="utf-8-sig"), delimiter="\t"
            ):
                if (
                    row["Project ID"] == "TCGA-LUAD"
                    and row["Sample Type"] in GROUPS
                    and row["File Name"].endswith("rna_seq.augmented_star_gene_counts.tsv")
                ):
                    fid = row["File ID"]
                    if fid in metadata and metadata[fid] != row:
                        raise ValueError("Conflicting metadata for File ID")
                    metadata[fid] = row
        with z.open("gdc_manifest.txt") as stream:
            manifest = {
                r["id"]: r
                for r in csv.DictReader(
                    io.TextIOWrapper(stream, encoding="utf-8-sig"), delimiter="\t"
                )
            }
        seen = set()
        with z.open("tcga_data.tar.gz") as compressed, tarfile.open(
            fileobj=compressed, mode="r|gz"
        ) as tar:
            for member in tar:
                if not member.isfile() or not member.name.endswith(
                    "rna_seq.augmented_star_gene_counts.tsv"
                ):
                    continue
                counts["expression_files_in_archive"] += 1
                parts = PurePosixPath(member.name).parts
                fid = parts[-2]
                if fid not in metadata:
                    counts["excluded_expression_files"] += 1
                    continue
                if fid in seen:
                    raise ValueError("Duplicate expression File ID in archive")
                seen.add(fid)
                row = metadata[fid]
                if parts[-1] != row["File Name"]:
                    raise ValueError("Sample-sheet filename mismatch")
                raw = tar.extractfile(member).read()
                if verify:
                    t = time.perf_counter()
                    m = manifest[fid]
                    if int(m["size"]) != len(raw) or hashlib.md5(raw).hexdigest() != m["md5"]:
                        raise ValueError(f"Manifest integrity failure: {fid}")
                    hash_s += time.perf_counter() - t
                    counts["md5_verified"] += 1
                t = time.perf_counter()
                tpm = get_tpm(raw, gene, method)
                parse_s += time.perf_counter() - t
                rows.append(
                    dict(
                        file_id=fid,
                        file_name=row["File Name"],
                        case_id=row["Case ID"],
                        sample_id=row["Sample ID"],
                        sample_type=row["Sample Type"],
                        gene=gene,
                        tpm=tpm,
                        log2_tpm_plus_1=math.log2(tpm + 1),
                    )
                )
                counts["eligible_files"] += 1
                if limit and len(rows) >= limit:
                    break
    counts.update(
        metadata_eligible_files=len(metadata),
        metadata_files_absent_from_processed_archive=len(set(metadata) - seen),
        parse_gene_seconds=parse_s,
        md5_seconds=hash_s,
    )
    if not rows:
        raise ValueError("No eligible gene expression records")
    return sorted(rows, key=lambda r: r["file_id"]), counts


def write_rows(path, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def run(args):
    started = time.perf_counter()
    rows, info = extract(
        args.archive, args.gene, args.method, args.limit, not args.skip_integrity_check
    )
    write_rows(args.output, rows)
    info.update(
        project=4,
        method=args.method,
        gene=args.gene,
        groups={group: sum(r["sample_type"] == group for r in rows) for group in sorted(GROUPS)},
        integrity_checked=not args.skip_integrity_check,
    )
    hint = [r["tpm"] for r in rows if r["file_id"] == "e6e83fd4-2f80-44e1-8dce-588a37308dbb"]
    if args.gene == "NKX2-1" and hint:
        if hint != [53.1718]:
            raise AssertionError("Problem-statement example does not match")
        info["statement_example_tpm"] = hint[0]
    json_write(Path(args.output).with_suffix(".summary.json"), info)
    emit_metrics(args, info, started)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--archive", default=str(ROOT / "data/tcga_source.zip"))
    p.add_argument("--gene", default="NKX2-1")
    p.add_argument("--method", choices=METHODS, default="bytes")
    p.add_argument("--output", default=str(ROOT / "results/tcga/NKX2-1.tsv"))
    p.add_argument(
        "--skip-integrity-check",
        action="store_true",
        help="Skip per-expression-file MD5; use same setting for all benchmarks",
    )
    p.add_argument(
        "--limit", type=int, default=0, help="Limit eligible expression files, not archive entries"
    )
    p.add_argument("--metrics")
    args = p.parse_args()
    if args.limit < 0 or not args.gene or "\t" in args.gene or "\n" in args.gene:
        p.error("Invalid gene or limit")
    run(args)
