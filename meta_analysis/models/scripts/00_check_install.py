"""Record the software and database versions used, and verify database provenance.

Run:  python3 scripts/00_check_install.py   (from ../models)
Writes outputs/environment.txt
"""
import hashlib
import os
import platform
import sys

import numpy
import pandas
import scipy
import matplotlib
import phreeqpython

sys.path.insert(0, os.path.dirname(__file__))
import iphreeqc  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DBDIR = os.path.join(ROOT, "databases_downloaded")
OUT = os.path.join(ROOT, "outputs", "environment.txt")

# md5 of the files as downloaded on 2026-10-07 from
# https://raw.githubusercontent.com/usgs-coupled/phreeqc3/master/database/<name>
EXPECTED_MD5 = {
    "phreeqc.dat": "07ba2ffc13d8e52eb35d4f992955e877",
    "wateq4f.dat": "d5d847569dcd20a19d4ae078910c7c60",
    "minteq.dat": "75e046ac2fe4f5a9afe2230e80cb1248",
    "minteq.v4.dat": "963c33b858001ff7be7121e9f38acd93",
    "llnl.dat": "b11e0d44aa3635395125420e7313a466",
    "Tipping_Hurley.dat": "9b2a682a37fe5ddcf08edcda2859630e",
}


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        h.update(f.read())
    return h.hexdigest()


lines = []
lines.append(f"python {sys.version.split()[0]} on {platform.platform()}")
lines.append(f"numpy {numpy.__version__}, scipy {scipy.__version__}, pandas {pandas.__version__}, matplotlib {matplotlib.__version__}")
lines.append(f"phreeqpython {phreeqpython.__version__ if hasattr(phreeqpython, '__version__') else '1.6.2 (pip show)'}")
lines.append(f"IPhreeqc library bundled with phreeqpython: {iphreeqc.version()}  ({iphreeqc._LIB_PATH})")
try:
    import importlib.metadata as _md
    import phreeqc as _phq
    lines.append(f"PyPI package phreeqc {_md.version('phreeqc')}: IPhreeqc {_phq.Phreeqc.GetVersionString()} (primary engine)")
except Exception as _e:  # pragma: no cover
    lines.append(f"PyPI package phreeqc not available: {_e}")
lines.append("")
lines.append("Databases (source: https://github.com/usgs-coupled/phreeqc3/tree/master/database, fetched 2026-10-07)")
for name, exp in EXPECTED_MD5.items():
    p = os.path.join(DBDIR, name)
    got = md5(p) if os.path.exists(p) else "MISSING"
    ok = "OK" if got == exp else "MISMATCH"
    # check that IPhreeqc 3.7.3 can load it
    try:
        iphreeqc.IPhreeqc(database_path=p)
        load = "loads in IPhreeqc"
    except Exception as e:  # pragma: no cover
        load = "LOAD FAILED: " + str(e).splitlines()[0]
    lines.append(f"  {name:20s} md5 {got} {ok}; {load}")
bundled = os.path.join(os.path.dirname(phreeqpython.__file__), "database")
lines.append("")
lines.append(f"Databases bundled inside phreeqpython ({bundled}): " + ", ".join(sorted(os.listdir(bundled))))
txt = "\n".join(lines)
print(txt)
os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w") as f:
    f.write(txt + "\n")
