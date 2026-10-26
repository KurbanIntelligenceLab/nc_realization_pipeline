"""Reproduces every estimate in the paper from the released derived data (cohorts and
matches in data_derived/) and compares it with the released results. Runs in seconds and
needs no download."""
from __future__ import annotations

import json
import math
import os
import shutil

import pytest

from ncr import config as C
from ncr.analysis import step_estimate
from ncr.util import read_json

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DD = os.path.join(ROOT, "data_derived")
pytestmark = pytest.mark.skipif(not os.path.exists(os.path.join(DD, "results", "summary.json")),
                                reason="released derived data not present")


def _compare(a, b, path, out):
    if isinstance(a, dict) and isinstance(b, dict):
        for k in a:
            if k in b and k not in ("versions", "figure_checks"):
                _compare(a[k], b[k], f"{path}/{k}", out)
    elif isinstance(a, list) and isinstance(b, list) and len(a) == len(b):
        for i, (x, y) in enumerate(zip(a, b)):
            _compare(x, y, f"{path}[{i}]", out)
    elif isinstance(a, (int, float)) and isinstance(b, (int, float)) and not isinstance(a, bool):
        out[0] += 1
        if not (abs(a - b) <= 1e-9 * max(1.0, abs(b)) or (math.isnan(a) and math.isnan(b))):
            out[1].append(f"{path}: {a} != {b}")


def test_released_estimates_are_reproduced(tmp_path):
    work = tmp_path / "work"
    for sub in ("cohorts", "matches", "checks"):
        shutil.copytree(os.path.join(DD, sub), work / sub)
    cfg_src = read_json(os.path.join(ROOT, "config", "study.json"))
    cfg_src["workdir"] = str(work)
    p = tmp_path / "study.json"
    p.write_text(json.dumps(cfg_src))
    new = step_estimate(C.load(str(p)))
    ref = read_json(os.path.join(DD, "results", "summary.json"))
    out = [0, []]
    for sec in ("base", "threshold", "capture_recapture"):
        _compare(new[sec], ref[sec], sec, out)
    assert out[0] > 500 and not out[1], "\n".join(out[1][:20])
