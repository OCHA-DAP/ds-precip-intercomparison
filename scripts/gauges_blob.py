"""Keep the raw rain-gauge downloads (data/gauges/) on the dev blob, verified by sha256.

data/ is git-ignored, so the blob is the shared copy of the gauge files and their
provenance (manifests, download logs, notification emails, scouting notes). `upload`
parks them with a sha256 manifest; `download` restores them and fails on any missing or
corrupt file. See docs/gauges/README.md for what each source is and how it was obtained.

    python scripts/gauges_blob.py upload     # data/gauges/ -> dev projects/<prefix>/raw/gauges/
    python scripts/gauges_blob.py download   # and back, sha256-verified

Neither direction deletes: `download` writes the manifest's files and leaves anything else
in data/gauges/ alone, and `upload` sends whatever is present locally (unchanged files are
skipped). To retire a file, delete it locally and on the blob, then upload.

TAHMO is NOT uploaded: its research licence says "Research use does not allow for sharing
of the data to any other person or organization", and the team blob is shared. It lives on
one machine; re-download it from the TAHMO portal and check it against
docs/gauges/tahmo.sha256 (ADR 0002).
"""

import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import ocha_stratus as stratus  # noqa: E402

from src.constants import PROJECT_PREFIX  # noqa: E402

ROOT = Path(__file__).resolve().parents[1] / "data" / "gauges"
PREFIX = f"{PROJECT_PREFIX}/raw/gauges"
MANIFEST = "manifest.json"
# top-level folders that must never leave this machine, with the reason
LOCAL_ONLY = {"tahmo": "TAHMO research licence: no sharing with any other person or organization (ADR 0002)"}
# empty hand-off markers between Claude sessions, not data
MARKERS = {"DONE", "DONE_round2"}


def _sha(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def _files() -> list[Path]:
    out = []
    for p in sorted(ROOT.rglob("*")):
        if not p.is_file() or p.name in MARKERS or p.name.startswith("."):
            continue
        if p.relative_to(ROOT).parts[0] in LOCAL_ONLY:
            continue
        out.append(p)
    return out


def upload() -> None:
    if not ROOT.is_dir():
        raise SystemExit(f"{ROOT} does not exist: nothing to upload")
    for name, why in LOCAL_ONLY.items():
        if (ROOT / name).exists():
            print(f"[gauges] not uploading {name}/: {why}")
    files = _files()
    if not files:
        raise SystemExit(f"no uploadable files under {ROOT}")
    cc = stratus.get_container_client("projects", stage="dev", write=True)
    manifest = {"created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                "local_only": LOCAL_ONLY, "files": {}}
    try:  # unchanged files (same sha256 as the current manifest) are not re-uploaded
        old = json.loads(cc.download_blob(f"{PREFIX}/{MANIFEST}").readall())["files"]
    except Exception:  # noqa: BLE001 - first upload: no manifest yet
        old = {}
    sent = 0
    for p in files:
        rel = p.relative_to(ROOT).as_posix()
        b = p.read_bytes()
        manifest["files"][rel] = {"sha256": _sha(b), "bytes": len(b)}
        if old.get(rel, {}).get("sha256") == manifest["files"][rel]["sha256"]:
            continue
        cc.upload_blob(f"{PREFIX}/{rel}", b, overwrite=True, max_concurrency=4)
        sent += 1
    # manifest last: an interrupted upload leaves the previous manifest, which then fails
    # verification on download instead of silently serving a mixed set
    cc.upload_blob(f"{PREFIX}/{MANIFEST}", json.dumps(manifest, indent=1).encode(), overwrite=True)
    print(f"[gauges] uploaded {sent} changed of {len(files)} files + manifest to dev projects/{PREFIX}")


def download() -> None:
    cc = stratus.get_container_client("projects", stage="dev")
    manifest = json.loads(cc.download_blob(f"{PREFIX}/{MANIFEST}").readall())
    bad = []
    for rel, info in manifest["files"].items():
        b = cc.download_blob(f"{PREFIX}/{rel}").readall()
        if _sha(b) != info["sha256"]:
            bad.append(rel)
            continue
        dest = ROOT / rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(b)
    if bad:
        raise SystemExit(f"[gauges] sha256 mismatch, not written: {bad}")
    print(f"[gauges] restored {len(manifest['files'])} files to {ROOT} (manifest {manifest['created']})")
    for name, why in manifest.get("local_only", {}).items():
        print(f"[gauges] {name}/ is not on blob ({why}); see docs/gauges/README.md to re-acquire")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("action", choices=["upload", "download"])
    {"upload": upload, "download": download}[ap.parse_args().action]()
