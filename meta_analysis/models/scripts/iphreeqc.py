"""Thin ctypes wrapper around the IPhreeqc shared library shipped with phreeqpython.

phreeqpython 1.6.2 (PyPI) bundles IPhreeqc 3.7.3-15968 as viphreeqc.so. Its own Python
wrapper does not expose the output, warning or selected-output strings, so this module
calls the library directly. Nothing here changes PHREEQC's numerics.
"""
import ctypes
import io
import os
import csv

import phreeqpython

_LIB_PATH = os.path.join(os.path.dirname(phreeqpython.__file__), "lib", "viphreeqc.so")
_lib = ctypes.cdll.LoadLibrary(_LIB_PATH)

_sig = {
    "CreateIPhreeqc": ([], ctypes.c_int),
    "DestroyIPhreeqc": ([ctypes.c_int], ctypes.c_int),
    "LoadDatabase": ([ctypes.c_int, ctypes.c_char_p], ctypes.c_int),
    "LoadDatabaseString": ([ctypes.c_int, ctypes.c_char_p], ctypes.c_int),
    "RunString": ([ctypes.c_int, ctypes.c_char_p], ctypes.c_int),
    "GetErrorString": ([ctypes.c_int], ctypes.c_char_p),
    "GetWarningString": ([ctypes.c_int], ctypes.c_char_p),
    "SetOutputStringOn": ([ctypes.c_int, ctypes.c_int], ctypes.c_int),
    "GetOutputString": ([ctypes.c_int], ctypes.c_char_p),
    "SetSelectedOutputStringOn": ([ctypes.c_int, ctypes.c_int], ctypes.c_int),
    "GetSelectedOutputString": ([ctypes.c_int], ctypes.c_char_p),
    "GetVersionString": ([], ctypes.c_char_p),
}
for _name, (_args, _res) in _sig.items():
    _f = getattr(_lib, _name)
    _f.argtypes = _args
    _f.restype = _res


def version():
    return _lib.GetVersionString().decode()


class IPhreeqc:
    def __init__(self, database_path=None, database_string=None):
        self.id = _lib.CreateIPhreeqc()
        if database_path is not None:
            n = _lib.LoadDatabase(self.id, str(database_path).encode())
        else:
            n = _lib.LoadDatabaseString(self.id, database_string.encode())
        if n != 0:
            raise RuntimeError("Database load errors:\n" + self.errors())
        _lib.SetSelectedOutputStringOn(self.id, 1)

    def errors(self):
        s = _lib.GetErrorString(self.id)
        return s.decode(errors="replace") if s else ""

    def warnings(self):
        s = _lib.GetWarningString(self.id)
        return s.decode(errors="replace") if s else ""

    def output(self):
        s = _lib.GetOutputString(self.id)
        return s.decode(errors="replace") if s else ""

    def run(self, input_string, keep_output=False):
        _lib.SetOutputStringOn(self.id, 1 if keep_output else 0)
        n = _lib.RunString(self.id, input_string.encode())
        if n != 0:
            raise RuntimeError(f"PHREEQC reported {n} error(s):\n" + self.errors())
        return self.selected_output()

    def selected_output(self):
        s = _lib.GetSelectedOutputString(self.id)
        if not s:
            return []
        txt = s.decode(errors="replace")
        rows = list(csv.reader(io.StringIO(txt), delimiter="\t"))
        if not rows:
            return []
        header = [h.strip() for h in rows[0]]
        out = []
        for r in rows[1:]:
            if not any(x.strip() for x in r):
                continue
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

    def __del__(self):
        try:
            _lib.DestroyIPhreeqc(self.id)
        except Exception:
            pass
