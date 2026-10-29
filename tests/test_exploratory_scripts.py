"""The three exploratory scripts run on the released data layout and reproduce the released results.

Each script reads data_derived/ as released (cohorts/, matches/ and config/study.json). A script that
expects any other layout fails here, as it would for a reader rerunning the paper."""
import json
import os
import shutil

import pytest

from scripts_path import *  # noqa: F401,F403  (adds scripts/ to the path)
import audit_estimates
import literature_sensitivity
import novelty_checks

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DD = os.path.join(ROOT, "data_derived")
RECORDS = os.path.join(ROOT, "data", "experiments", "experimental_records.csv.gz")
pytestmark = pytest.mark.skipif(
    not all(os.path.exists(os.path.join(DD, d, "oqmd.csv")) for d in ("cohorts", "matches")),
    reason="released derived data not present")


def _layout(tmp_path):
    d = tmp_path / "derived"
    for sub in ("cohorts", "matches"):
        shutil.copytree(os.path.join(DD, sub), d / sub)
    (d / "checks").mkdir()
    return d


def _leaves(x, p=""):
    if isinstance(x, dict):
        for k, v in x.items():
            yield from _leaves(v, f"{p}/{k}")
    elif isinstance(x, list):
        for i, v in enumerate(x):
            yield from _leaves(v, f"{p}[{i}]")
    else:
        yield p, x


def _assert_reproduces(d, name):
    new = dict(_leaves(json.load(open(d / "checks" / f"{name}.json"))))
    pub = dict(_leaves(json.load(open(os.path.join(DD, "checks", f"{name}.json")))))
    assert set(new) == set(pub)
    for k, v in pub.items():
        if isinstance(v, float):
            assert abs(new[k] - v) <= 1e-9, k
        else:
            assert new[k] == v, k


def test_audit_estimates_reproduces(tmp_path):
    d = _layout(tmp_path)
    audit_estimates.main(["--derived", str(d), "--audit", os.path.join(ROOT, "data", "audit", "prior_record_audit.csv")])
    _assert_reproduces(d, "audit_estimates")


def test_literature_sensitivity_reproduces(tmp_path):
    d = _layout(tmp_path)
    literature_sensitivity.main(["--derived", str(d), "--audit", os.path.join(ROOT, "data", "audit", "literature_audit.csv")])
    _assert_reproduces(d, "literature_sensitivity")


@pytest.mark.skipif(not os.path.exists(RECORDS), reason="experimental records not present")
def test_novelty_checks_reproduces(tmp_path):
    d = _layout(tmp_path)
    novelty_checks.main(["--derived", str(d), "--records", RECORDS])
    _assert_reproduces(d, "novelty_checks")
