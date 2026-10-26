"""Checks on the released inputs and derived data: the inputs match the checksums recorded in
the study, and the released results match the numbers reported in the paper."""
from __future__ import annotations

import gzip
import hashlib
import importlib.util
import json
import os

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _sha256_gz(path):
    h = hashlib.sha256()
    with gzip.open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _checksums():
    out = []
    with open(os.path.join(ROOT, "data", "CHECKSUMS.txt"), encoding="utf-8") as fh:
        for line in fh:
            if line.strip() and not line.startswith("#"):
                digest, rel = line.split()
                out.append((digest, os.path.join(ROOT, "data", rel + ".gz")))
    return out


@pytest.mark.parametrize("digest,path", _checksums())
def test_inputs_match_the_study_checksums(digest, path):
    if not os.path.exists(path):
        pytest.skip(f"{os.path.relpath(path, ROOT)} not present")
    assert _sha256_gz(path) == digest


def test_released_results_match_the_paper():
    work = os.path.join(ROOT, "data_derived")
    if not os.path.exists(os.path.join(work, "results", "robustness.json")):
        pytest.skip("released derived data not present")
    spec = importlib.util.spec_from_file_location("compare_results", os.path.join(ROOT, "scripts", "compare_results.py"))
    cr = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cr)
    with open(os.path.join(ROOT, "tests", "paper_numbers.json"), encoding="utf-8") as fh:
        ref = json.load(fh)
    assert cr.diff(cr.collect(work), ref) == []
