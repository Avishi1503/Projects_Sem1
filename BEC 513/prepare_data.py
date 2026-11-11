"""Prepare the input files from the supplied local archives."""
import argparse
import hashlib
import json
import shutil
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--downloads", type=Path, default=Path("D:/Downloads"))
    args = p.parse_args()
    data = ROOT / "data"
    data.mkdir(exist_ok=True)
    source = args.downloads / "27087241.zip"
    if source.exists():
        with zipfile.ZipFile(source) as z:
            for name in [
                "shuf.a.bed.gz",
                "shuf.b.bed.gz",
                "query.bed.gz",
                "reference.hist",
                "selection.tsv",
            ]:
                target = data / name
                if target.exists() and target.stat().st_size == z.getinfo(name).file_size:
                    continue
                with z.open(name) as f, target.open("wb") as g:
                    shutil.copyfileobj(f, g, 8 * 1024 * 1024)
    for name, target_name in [
        ("mapped.bed.gz", "mapped.bed.gz"),
        ("27129882.zip", "tcga_source.zip"),
    ]:
        source = args.downloads / name
        target = data / target_name
        if source.exists() and (
            not target.exists() or target.stat().st_size != source.stat().st_size
        ):
            shutil.copy2(source, target)
    required = [
        "shuf.a.bed.gz",
        "shuf.b.bed.gz",
        "query.bed.gz",
        "reference.hist",
        "selection.tsv",
        "mapped.bed.gz",
        "tcga_source.zip",
    ]
    for name in required:
        if not (data / name).exists():
            raise FileNotFoundError(f"Missing {data/name}; supply the original local file")
    records = [
        dict(path=name, bytes=(data / name).stat().st_size, sha256=digest(data / name))
        for name in required
    ]
    (data / "prepared_inputs.json").write_text(json.dumps(records, indent=2))
    print("Prepared and hashed all local inputs. No download performed.")


if __name__ == "__main__":
    main()
