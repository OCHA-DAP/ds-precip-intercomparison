"""Move the generated site outputs (report data/figures, explorer data) to and from blob.

Generated JSON/PNGs never enter git: the analysis runs locally (or on
Databricks), `upload` parks them on the dev blob with a sha256 manifest, and the
Pages deploy workflow runs `download`, which verifies every file against the
manifest — a missing or corrupt file is a red deploy, not a gutted site.

    python scripts/site_blob.py upload
    python scripts/site_blob.py download
"""

import argparse
import hashlib
import json
import mimetypes
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import ocha_stratus as stratus  # noqa: E402
from azure.storage.blob import ContentSettings  # noqa: E402

from src.constants import PROJECT_PREFIX  # noqa: E402

ROOT = Path(__file__).resolve().parents[1] / "site"
SUBDIRS = ("report/data", "report/fig", "explorer/data")
PREFIX = f"{PROJECT_PREFIX}/site"
MANIFEST = "manifest.json"


def _sha(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def upload() -> None:
    files = sorted(p for sub in SUBDIRS for p in (ROOT / sub).rglob("*") if p.is_file())
    if not files:
        raise SystemExit("nothing to upload (run analyze.py and figures.py first)")
    cc = stratus.get_container_client("projects", stage="dev", write=True)
    manifest = {"created": datetime.now(timezone.utc).isoformat(timespec="seconds"), "files": {}}
    # unchanged files (same sha256 as the current manifest) are not re-uploaded
    try:
        old = json.loads(cc.download_blob(f"{PREFIX}/{MANIFEST}").readall())["files"]
    except Exception:  # noqa: BLE001 - first upload
        old = {}
    sent = 0
    for p in files:
        rel = p.relative_to(ROOT).as_posix()
        b = p.read_bytes()
        sha = _sha(b)
        manifest["files"][rel] = sha
        if old.get(rel) == sha:
            continue
        ctype = mimetypes.guess_type(p.name)[0] or "application/octet-stream"
        cc.upload_blob(f"{PREFIX}/{rel}", b, overwrite=True, content_settings=ContentSettings(content_type=ctype),
                       max_concurrency=4)
        sent += 1
    # manifest last: a partial upload leaves the previous manifest, which then
    # fails verification on download instead of serving a mixed set
    cc.upload_blob(f"{PREFIX}/{MANIFEST}", json.dumps(manifest, indent=1).encode(), overwrite=True,
                   content_settings=ContentSettings(content_type="application/json"))
    print(f"uploaded {sent} changed of {len(files)} files + manifest to projects/{PREFIX}")


def download() -> None:
    cc = stratus.get_container_client("projects", stage="dev")
    manifest = json.loads(cc.download_blob(f"{PREFIX}/{MANIFEST}").readall())
    for sub in SUBDIRS:
        shutil.rmtree(ROOT / sub, ignore_errors=True)
    bad = []
    for rel, sha in manifest["files"].items():
        b = cc.download_blob(f"{PREFIX}/{rel}").readall()
        if _sha(b) != sha:
            bad.append(rel)
            continue
        out = ROOT / rel
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_bytes(b)
    if bad:
        raise SystemExit(f"{len(bad)} files failed sha256 verification, e.g. {bad[:3]}")
    print(f"downloaded {len(manifest['files'])} verified files (manifest {manifest['created']})")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("action", choices=["upload", "download"])
    upload() if ap.parse_args().action == "upload" else download()
