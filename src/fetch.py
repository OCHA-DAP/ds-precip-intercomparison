"""HTTP download helpers (retries, gzip, NASA Earthdata auth)."""

import gzip
import os
import shutil
import time
from pathlib import Path

import requests

URS_HOST = "urs.earthdata.nasa.gov"


class AuthError(RuntimeError):
    """Login/approval page instead of data: do not retry."""


class EarthdataSession(requests.Session):
    """Keeps basic auth across the GES DISC -> URS -> GES DISC redirect chain
    (requests drops it on cross-host redirects by default)."""

    def __init__(self, username: str, password: str):
        super().__init__()
        self.auth = (username, password)

    def rebuild_auth(self, prepared_request, response):
        headers = prepared_request.headers
        url = prepared_request.url
        if "Authorization" in headers:
            orig = requests.utils.urlparse(response.request.url).hostname
            redirect = requests.utils.urlparse(url).hostname
            if orig != redirect and redirect != URS_HOST and orig != URS_HOST:
                del headers["Authorization"]


def earthdata_session() -> requests.Session:
    user = os.environ.get("IMERG_USERNAME") or os.environ.get("EARTHDATA_USERNAME")
    pw = os.environ.get("IMERG_PASSWORD") or os.environ.get("EARTHDATA_PASSWORD")
    if user and pw:
        return EarthdataSession(user, pw)
    # Fall back to ~/.netrc (requests reads it for the URS host on redirect).
    return requests.Session()


def download(
    url: str,
    dest: Path,
    session: requests.Session | None = None,
    retries: int = 8,
    timeout: int = 300,
    gunzip: bool = False,
) -> Path:
    """Stream `url` to `dest` (atomic rename). gunzip=True stores the
    decompressed payload. Skips if dest already exists."""
    dest = Path(dest)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    sess = session or requests
    tmp = dest.with_suffix(dest.suffix + ".part")
    raw = tmp.with_suffix(".gz") if gunzip else tmp
    raw.unlink(missing_ok=True)
    last = None
    for attempt in range(retries):
        try:
            # Resume a broken transfer only if the server answers 206 to a Range
            # request; a 200 means it ignored Range, so start over (never append
            # a full body onto a partial file).
            have = raw.stat().st_size if raw.exists() else 0
            headers = {"Range": f"bytes={have}-"} if have else {}
            with sess.get(url, stream=True, timeout=timeout, headers=headers) as r:
                r.raise_for_status()
                if "text/html" in r.headers.get("Content-Type", ""):
                    raise AuthError(f"got HTML instead of data (auth/approval page?) for {url}")
                mode = "ab" if (have and r.status_code == 206) else "wb"
                expected = r.headers.get("Content-Length")
                with open(raw, mode) as f:
                    for chunk in r.iter_content(chunk_size=1 << 20):
                        f.write(chunk)
            if expected is not None and r.status_code in (200, 206):
                total = int(expected) + (have if mode == "ab" else 0)
                if raw.stat().st_size != total:
                    raise IOError(f"short read {raw.stat().st_size} != {total}")
            if gunzip:
                with gzip.open(raw, "rb") as fin, open(tmp if raw != tmp else dest, "wb") as fout:
                    shutil.copyfileobj(fin, fout)
                raw.unlink()
            tmp.rename(dest)
            return dest
        except AuthError:
            raw.unlink(missing_ok=True)
            raise
        except Exception as exc:  # noqa: BLE001 - retry anything transient
            last = exc
            code = exc.response.status_code if isinstance(exc, requests.HTTPError) and exc.response is not None else None
            # Never retry auth failures: repeated bad logins lock the Earthdata
            # account (shared with the daily IMERG pipeline) for 10 minutes.
            if code in (401, 403, 404):
                raw.unlink(missing_ok=True)
                raise
            if isinstance(exc, (OSError, gzip.BadGzipFile)) and gunzip and not isinstance(exc, requests.RequestException):
                raw.unlink(missing_ok=True)  # corrupt archive: restart from scratch
            time.sleep(min(60, 5 * 2**attempt))
    raw.unlink(missing_ok=True)
    raise RuntimeError(f"download failed after {retries} attempts: {url}: {last}")


def exists(url: str, session: requests.Session | None = None) -> bool:
    sess = session or requests
    try:
        r = sess.head(url, timeout=60, allow_redirects=True)
        return r.status_code == 200
    except requests.RequestException:
        return False
