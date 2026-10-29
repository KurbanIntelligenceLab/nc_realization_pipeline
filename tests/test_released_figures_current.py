"""The Source Data stored in data_derived/figures/ is current: regenerating it from the released
figure data and audit gives the same workbooks, cell for cell, formulas included. A figure revision
that is not re-exported fails here."""
import json
import os

import pytest

from ncr import source_data

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DD = os.path.join(ROOT, "data_derived")
FD = os.path.join(DD, "results", "figure_data.json")
SD = os.path.join(DD, "figures", "source_data")
pytestmark = pytest.mark.skipif(not (os.path.exists(FD) and os.path.isdir(SD)), reason="released derived data not present")


def _cells(path):
    import openpyxl
    wb = openpyxl.load_workbook(path, data_only=False)
    return {(ws.title, c.coordinate): c.value for ws in wb.worksheets for row in ws.iter_rows() for c in row
            if c.value is not None}


def test_stored_source_data_matches_regenerated(tmp_path):
    figs = json.load(open(FD))
    audit = source_data.read_audit(os.path.join(ROOT, "data", "audit", "literature_audit.csv"))
    written = source_data.write_all(figs, str(tmp_path), audit)
    assert written
    for name in sorted(os.listdir(tmp_path)):
        if not name.endswith(".xlsx"):
            continue
        new, old = _cells(tmp_path / name), _cells(os.path.join(SD, name))
        assert set(new) == set(old), name
        for k, v in new.items():
            w = old[k]
            if isinstance(v, float) and isinstance(w, (int, float)):
                assert abs(v - w) <= 1e-9 * max(1.0, abs(v)), (name, k)
            else:
                assert v == w, (name, k)
