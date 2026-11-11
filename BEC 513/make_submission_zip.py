"""Package code, reports and compact results; omit supplied bulk data and BEDs."""
import hashlib
import json
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent
destination = ROOT / "Python_Bioinformatics_Submission.zip"
data_allowed = {"selection.tsv", "reference.hist", "sources.json", "prepared_inputs.json"}
included = []
with zipfile.ZipFile(destination, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as z:
    for path in sorted(ROOT.rglob("*")):
        if not path.is_file():
            continue
        rel = path.relative_to(ROOT)
        if any(x in rel.parts for x in ("__pycache__", ".git", ".venv", "tmp")):
            continue
        if path.suffix.lower() in (".zip", ".gz", ".bed", ".npy", ".nbc", ".nbi", ".pyc"):
            continue
        if rel.parts[0] == "data" and (len(rel.parts) != 2 or path.name not in data_allowed):
            continue
        z.write(path, arcname="Python_Bioinformatics_Projects/" + rel.as_posix())
        included.append(dict(path=rel.as_posix(), bytes=path.stat().st_size))
with zipfile.ZipFile(destination) as z:
    if z.testzip() is not None:
        raise RuntimeError("ZIP integrity check failed")
print(
    json.dumps(
        dict(zip=str(destination), bytes=destination.stat().st_size, files=len(included)), indent=2
    )
)
