"""Common interface to two independent builds of the USGS IPhreeqc engine.

* "386": IPhreeqc 3.8.6-17100, from the PyPI package `phreeqc` 1.1.1
  (pip install phreeqc --break-system-packages). Primary engine for all results.
* "373": IPhreeqc 3.7.3-15968, bundled with the PyPI package `phreeqpython` 1.6.2
  (pip install phreeqpython --break-system-packages), called through ctypes in iphreeqc.py.
  Used only as an engine cross-check. It cannot read the 2024 phreeqc.dat (new keyword
  GAS_BINARY_PARAMETERS) but reads minteq.v4.dat and wateq4f.dat.

Both return SELECTED_OUTPUT as a list of dicts (one dict per selected-output row).
"""
import csv
import io
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DBDIR = os.path.join(ROOT, "databases_downloaded")


def db_path(name):
    return os.path.join(DBDIR, name)


def _is_header(r):
    """A header row has no numeric fields (PHREEQC may emit a stale header from an earlier run)."""
    vals = [x.strip() for x in r if x.strip()]
    if not vals:
        return False
    for v in vals:
        try:
            float(v)
            return False
        except ValueError:
            pass
    return True


def _parse_selected(txt):
    rows = list(csv.reader(io.StringIO(txt), delimiter="\t"))
    if not rows:
        return []
    header = None
    out = []
    for r in rows:
        if not any(x.strip() for x in r):
            continue
        if _is_header(r):
            header = [h.strip() for h in r]   # always use the most recent header
            continue
        if header is None:
            raise RuntimeError("selected output without header")
        if len(r) != len(header):
            raise RuntimeError(f"selected output row has {len(r)} fields, header {len(header)}")
        rec = {}
        for h, v in zip(header, r):
            if h == "":
                continue
            v = v.strip()
            try:
                rec[h] = float(v)
            except ValueError:
                rec[h] = v
        out.append(rec)
    return out


class Engine:
    def __init__(self, database, which="386"):
        self.which = which
        self.database = database
        if which == "386":
            import phreeqc
            self._p = phreeqc.Phreeqc()
            n = self._p.LoadDatabase(database)
            if n != 0:
                raise RuntimeError("database load errors:\n" + self._p.GetErrorString())
            self._p.SetSelectedOutputStringOn(True)
            self.version = phreeqc.Phreeqc.GetVersionString()
        elif which == "373":
            import iphreeqc
            self._p = iphreeqc.IPhreeqc(database_path=database)
            self.version = iphreeqc.version()
        else:
            raise ValueError(which)

    def run(self, text, keep_output=False):
        if self.which == "386":
            self._p.SetOutputStringOn(bool(keep_output))
            n = self._p.RunString(text)
            if n != 0:
                raise RuntimeError(f"PHREEQC {self.version}: {n} error(s)\n" + self._p.GetErrorString())
            self.last_output = self._p.GetOutputString() if keep_output else ""
            self.last_warnings = self._p.GetWarningString()
            return _parse_selected(self._p.GetSelectedOutputString())
        else:
            import iphreeqc
            self._p.run(text, keep_output=keep_output)
            self.last_output = self._p.output() if keep_output else ""
            self.last_warnings = self._p.warnings()
            raw = iphreeqc._lib.GetSelectedOutputString(self._p.id)
            return _parse_selected(raw.decode(errors="replace") if raw else "")
