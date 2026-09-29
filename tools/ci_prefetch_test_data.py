#!/usr/bin/env python3
"""Pre-seed the astropy download cache with the remote data files used by the
unit-test suite and the tutorial-notebook tests.

Those pull data files from the pycbc-config repo and the pycbc_data release
assets via ``pycbc.io.get_file`` / ``astropy.utils.data.download_file(...,
cache=True)``.  On a cold CI cache every job re-downloads them, and a
transient raw.githubusercontent.com / GitHub error then fails the job.

Running this script populates ``~/.astropy/cache/download`` directly (it only
needs the Python standard library, not pycbc or astropy).  That directory is
persisted with ``actions/cache`` and refreshed daily by the
``refresh-ci-caches`` workflow, so test runs get cache hits instead.

Keep ``URLS`` in step with the tests that download with ``cache=True``.
Downloading is best effort: ``get_file`` still retries at test time, and a
miss just means the test fetches the file itself.
"""

import hashlib
import os
import sys
import time
import urllib.parse
import urllib.request

PYCBC_CONFIG = (
    "https://github.com/gwastro/pycbc-config/raw/master/test_data_files/{}"
)
RELEASE = "https://github.com/gwastro/pycbc_data/releases/download/ci-data/{}"

URLS = []

# test/test_frame.py
URLS += [
    PYCBC_CONFIG.format(f"frametest{dtype}.gwf")
    for dtype in ("float32", "float64", "complex64", "complex128")
]

# test/test_skymax.py  (idx, jdx both range(4))
URLS += [
    PYCBC_CONFIG.format(f"skymaxtest_{kind}_{i}.hdf")
    for kind in ("stilde", "hplus", "hcross")
    for i in range(4)
]

# test/test_tmpltbank.py
URLS += [
    PYCBC_CONFIG.format(name)
    for name in (
        "ZERO_DET_high_P.txt.gz",
        "stockEvals.dat.gz",
        "stockEvecs.dat.gz",
        "stockChirps.dat.gz",
        "stockHexagonal.dat.gz",
        "stockAnstar3D.dat.gz",
    )
]

# test/test_live_coinc_compare.py
URLS += [PYCBC_CONFIG.format("H1L1-PTA_HISTOGRAM.hdf")]

# test/test_infmodel.py
URLS += [
    RELEASE.format(f"{ifo}-{ifo}1_LOSC_CLN_4_V1-1187007040-2048.gwf")
    for ifo in "HLV"
]

# test/test_bhspec.py
URLS += [
    RELEASE.format(f"{ifo}-{ifo}1_GWOSC_4KHZ_R1-1126259447-32.gwf")
    for ifo in "HL"
]

# PyCBC-Tutorials notebooks run in the "tutorial tests" workflow.  They pass
# get_file the original GWOSC / DCC URL, which (inside Actions) is fetched from
# the release asset -- so key the cache on the release URL.  The LOSC_CLN 2048s
# files inference_0 / inference_1 use are already listed above for
# test/test_infmodel.py.
URLS += [
    RELEASE.format(name)
    for name in (
        "H-H1_LOSC_4_V2-1126259446-32.gwf",   # examples/gw150914_audio
        "L-L1_LOSC_4_V2-1135136334-32.gwf",   # examples/gw151226_snr
        "H-H1_LOSC_4_V1-1167559920-32.gwf",   # examples/gw170104_look
        "L-L1_LOSC_4_V1-1167559920-32.gwf",   # examples/gw170104_look
        "H-H1_LOSC_4_V2-1128678884-32.gwf",   # tutorial/1_CatalogData
    )
]
URLS += [
    RELEASE.format(f"{ifo}-{ifo}1_GWOSC_4KHZ_R2-1242440920-4096.gwf")
    for ifo in "HLV"                          # tutorial/inference_8_BHRingdown
]

CACHE_DIR = os.path.join(
    os.path.expanduser("~"), ".astropy", "cache", "download", "url"
)


def _url_to_dirname(url):
    """Reproduce astropy.utils.data._url_to_dirname."""
    parts = list(urllib.parse.urlsplit(url))
    parts[1] = parts[1].lower()
    if parts[0].lower() in ("http", "https") and parts[1] and parts[2] == "":
        parts[2] = "/"
    normalised = urllib.parse.urlunsplit(parts)
    return hashlib.md5(normalised.encode("utf-8")).hexdigest()


def prefetch(url, attempts=5):
    entry = os.path.join(CACHE_DIR, _url_to_dirname(url))
    contents = os.path.join(entry, "contents")
    if os.path.exists(contents):
        return "cached"
    os.makedirs(entry, exist_ok=True)
    last_exc = None
    for i in range(1, attempts + 1):
        try:
            req = urllib.request.Request(
                url, headers={"User-Agent": "pycbc-ci-prefetch"}
            )
            with urllib.request.urlopen(req, timeout=30) as resp:
                data = resp.read()
            tmp = contents + ".tmp"
            with open(tmp, "wb") as fh:
                fh.write(data)
            os.replace(tmp, contents)
            with open(os.path.join(entry, "url"), "w") as fh:
                fh.write(url)
            return "ok"
        except Exception as exc:  # best effort
            last_exc = exc
            if i < attempts:
                time.sleep(min(5, 2 ** i))
    raise last_exc


def main():
    failures = []
    for url in URLS:
        try:
            status = prefetch(url)
            print(f"{status:6s} {url}")
        except Exception as exc:
            failures.append(url)
            print(f"FAIL   {url}  ({exc})")
    if failures:
        print(f"\n{len(failures)}/{len(URLS)} downloads failed")
        return 1
    print(f"\nall {len(URLS)} test-data files present in cache")
    return 0


if __name__ == "__main__":
    sys.exit(main())
