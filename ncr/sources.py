"""Dated experimental-structure sources derived from ICSD that are openly available, and the
merge of all sources into the single experimental-reports table read by `dbmatch`.

  OQMD-ICSD  ICSD-origin entries of an OQMD release with the publication year of their
             reference (in v1.8 the years end in 2012)
  MP-ICSD    Materials Project materials that provenance marks non-theoretical AND that carry
             ICSD identifiers and dated references; the date is the earliest non-default
             reference year (in the 2022-10-28 and 2025-09-25 builds the years end in 2016).
             Non-theoretical materials without any identifier or reference are NOT used: in
             the 2025 build 200,487 such records carry no evidence of an experimental match.
"""
from __future__ import annotations

import collections
import csv
import glob
import gzip
import json
import os
import re
import sys

from .util import log, write_csv, write_json

csv.field_size_limit(sys.maxsize)
OUT_FIELDS = ["exp_id", "database", "formula", "first_report_date", "experimental", "structure_file",
              "high_pressure", "thin_film_only"]
_YEAR = re.compile(r"year\s*=\s*[{\"]?\s*(\d{4})")


def convert_oqmd_icsd(in_tsv: str, out_csv: str) -> dict:
    """Input: tab-separated export with columns entry_id, composition_id, year."""
    from pymatgen.core import Composition
    rep, rows = collections.Counter(), []
    with open(in_tsv, newline="", encoding="utf-8", errors="replace") as fh:
        for r in csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE):
            rep["rows_read"] += 1
            y, comp = (r.get("year") or "").strip(), (r.get("composition_id") or "").strip()
            if not y.isdigit() or comp in ("", "NULL"):
                rep["excluded_no_year_or_formula"] += 1
                continue
            try:
                f = Composition(comp.replace(" ", "")).formula.replace(" ", "")
            except Exception:  # noqa: BLE001
                rep["excluded_unparseable_formula"] += 1
                continue
            rows.append({"exp_id": f"oqmdicsd-{r['entry_id']}", "database": "OQMD-ICSD", "formula": f,
                         "first_report_date": y, "experimental": 1, "structure_file": "", "high_pressure": "",
                         "thin_film_only": ""})
    years = [int(r["first_report_date"]) for r in rows]
    out = dict(rep, rows_written=len(rows), year_max=max(years) if years else None)
    write_csv(out_csv, rows, OUT_FIELDS)
    write_json(out_csv.replace(".csv", "_conversion_report.json"), out)
    log(f"convert-oqmd-icsd: {out}")
    return out


def convert_mp_icsd(build_dir: str, out_csv: str, thermo_type: str = "GGA_GGA+U") -> dict:
    from pymatgen.core import Composition
    docs = []
    for p in glob.glob(os.path.join(build_dir, "provenance", "**", "*.jsonl.gz"), recursive=True):
        with gzip.open(p, "rt", encoding="utf-8") as fh:
            for ln in fh:
                if ln.strip():
                    d = json.loads(ln)
                    docs.append((d.get("material_id"), bool(d.get("theoretical")), d.get("references") or [],
                                 (d.get("database_IDs") or {}).get("icsd") or []))
    default = collections.Counter()
    for _, th, refs, _ in docs:
        if th:
            default.update(set(refs))
    common = {r for r, c in default.items() if c >= 100}
    formula = {}
    for p in glob.glob(os.path.join(build_dir, "thermo", "**", "*.jsonl.gz"), recursive=True):
        with gzip.open(p, "rt", encoding="utf-8") as fh:
            for ln in fh:
                if not ln.strip():
                    continue
                d = json.loads(ln)
                mid = d.get("material_id")
                if mid in formula or not d.get("composition_reduced"):
                    continue
                formula[mid] = "".join(f"{el}{amt:g}" for el, amt in d["composition_reduced"].items())
    rep, rows = collections.Counter(), []
    for mid, th, refs, icsd in docs:
        if th:
            continue
        rep["non_theoretical"] += 1
        if not icsd:
            rep["excluded_no_icsd_identifier"] += 1
            continue
        ys = [int(m.group(1)) for r in refs if r not in common for m in [_YEAR.search(r)] if m]
        if not ys:
            rep["excluded_no_dated_reference"] += 1
            continue
        if mid not in formula:
            rep["excluded_no_formula"] += 1
            continue
        try:
            f = Composition(formula[mid]).formula.replace(" ", "")
        except Exception:  # noqa: BLE001
            rep["excluded_unparseable_formula"] += 1
            continue
        rows.append({"exp_id": f"mpicsd-{mid}", "database": "MP-ICSD", "formula": f, "first_report_date": str(min(ys)),
                     "experimental": 1, "structure_file": "", "high_pressure": "", "thin_film_only": ""})
    years = [int(r["first_report_date"]) for r in rows]
    out = dict(rep, rows_written=len(rows), default_reference_strings=len(common),
               year_max=max(years) if years else None)
    write_csv(out_csv, rows, OUT_FIELDS)
    write_json(out_csv.replace(".csv", "_conversion_report.json"), out)
    log(f"convert-mp-icsd: {out}")
    return out


def merge_experiments(inputs: list[str], out_csv: str) -> dict:
    rows, counts = [], collections.Counter()
    for p in inputs:
        with open(p, newline="", encoding="utf-8") as fh:
            for r in csv.DictReader(fh):
                rows.append(r)
                counts[r.get("database", "?")] += 1
    fields = []
    for r in rows:
        for k in r:
            if k not in fields:
                fields.append(k)
    write_csv(out_csv, rows, fields)
    out = {"inputs": inputs, "rows": len(rows), "by_database": dict(counts)}
    write_json(out_csv.replace(".csv", "_merge_report.json"), out)
    log(f"merge-experiments: {out}")
    return out
