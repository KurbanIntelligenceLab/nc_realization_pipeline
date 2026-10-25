"""COD (Crystallography Open Database): export from its public read-only MySQL server and
conversion to the experimental-reports table read by `dbmatch`.

Rules, fixed from the `data` table listing (535,176 rows when checked):
  formula       `formula` with the surrounding dashes removed, falling back to `calcformula`
  date          `year` is the publication year; year-only dates are set to 1 July (Table S1).
                COD's `date` column is the last modification date and is not used
  experimental  rows with method = 'theoretical' or status = 'retracted' are excluded
  high pressure cellpressure or diffrpressure above 1 GPa (1e6 kPa); an implementation choice
  disorder      flags contains 'has disorder'
CIF files are not exported here; `fetch-cod-cifs` downloads them only for matched entries.
"""
from __future__ import annotations

import csv
import collections
import os
import sys

from .util import log, write_csv, write_json

FIELDS = ["file", "formula", "calcformula", "nel", "sgNumber", "sg", "a", "b", "c", "alpha", "beta", "gamma",
          "vol", "Z", "year", "journal", "doi", "method", "status", "flags", "duplicateof", "onhold",
          "cellpressure", "diffrpressure", "celltemp", "mineral", "compoundsource"]
HIGH_PRESSURE_KPA = 1.0e6
csv.field_size_limit(sys.maxsize)


def export_cod(out_tsv: str, host: str = "sql.crystallography.net", batch: int = 5000) -> int:
    import pymysql
    import pymysql.cursors
    os.makedirs(os.path.dirname(os.path.abspath(out_tsv)), exist_ok=True)
    con = pymysql.connect(host=host, user="cod_reader", password="", database="cod", connect_timeout=30,
                          read_timeout=3600, cursorclass=pymysql.cursors.SSCursor)
    n = 0
    tmp = out_tsv + ".part"
    with con.cursor() as c, open(tmp, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", quoting=csv.QUOTE_MINIMAL)
        w.writerow(FIELDS)
        c.execute("SELECT " + ", ".join(f"`{f}`" for f in FIELDS) + " FROM data")
        while True:
            rows = c.fetchmany(batch)
            if not rows:
                break
            for r in rows:
                w.writerow(["" if v is None else str(v).replace("\t", " ").replace("\n", " ") for v in r])
            n += len(rows)
            if n % 50000 < batch:
                log(f"export-cod: {n} rows")
    con.close()
    os.replace(tmp, out_tsv)
    log(f"export-cod: {n} rows written to {out_tsv}")
    return n


def _clean_formula(s: str) -> str:
    return (s or "").strip().strip("-").strip()


def _f(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def convert_cod(in_tsv: str, out_csv: str) -> dict:
    from pymatgen.core import Composition
    import warnings
    warnings.filterwarnings("ignore")
    rep = collections.Counter()
    years = collections.Counter()
    rows = []
    with open(in_tsv, newline="", encoding="utf-8", errors="replace") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rep["rows_read"] += 1
            if r.get("method") == "theoretical":
                rep["excluded_theoretical"] += 1
                continue
            if r.get("status") == "retracted":
                rep["excluded_retracted"] += 1
                continue
            f = _clean_formula(r.get("formula")) or _clean_formula(r.get("calcformula"))
            if not f:
                rep["excluded_no_formula"] += 1
                continue
            try:
                comp = Composition(f.replace(" ", ""))
                if not comp.elements:
                    raise ValueError("empty")
            except Exception:  # noqa: BLE001
                rep["excluded_unparseable_formula"] += 1
                continue
            y = r.get("year", "")
            if not y or not y.isdigit() or int(y) < 1800:
                rep["excluded_no_year"] += 1
                continue
            years[int(y)] += 1
            els = {str(e) for e in comp.elements}
            if "C" in els and "H" in els:
                rep["contains_C_and_H"] += 1
            p = max(_f(r.get("cellpressure")) or 0.0, _f(r.get("diffrpressure")) or 0.0)
            rows.append({"exp_id": f"cod-{r['file']}", "database": "COD", "formula": comp.formula.replace(" ", ""),
                         "first_report_date": y, "experimental": 1, "structure_file": "",
                         "high_pressure": int(p > HIGH_PRESSURE_KPA), "thin_film_only": "",
                         "disordered": int("has disorder" in (r.get("flags") or "")),
                         "sg_number": r.get("sgNumber", ""), "doi": r.get("doi", ""),
                         "duplicateof": r.get("duplicateof", ""), "cod_file": r["file"]})
    rep["rows_written"] = len(rows)
    rep["high_pressure_rows"] = sum(x["high_pressure"] for x in rows)
    rep["disordered_rows"] = sum(x["disordered"] for x in rows)
    rep["rows_without_both_C_and_H"] = len(rows) - rep["contains_C_and_H"]
    out = dict(rep)
    out["year_min"], out["year_max"] = (min(years), max(years)) if years else (None, None)
    out["rows_by_year_2010_onward"] = {str(y): years[y] for y in sorted(years) if y >= 2010}
    fields = ["exp_id", "database", "formula", "first_report_date", "experimental", "structure_file", "high_pressure",
              "thin_film_only", "disordered", "sg_number", "doi", "duplicateof", "cod_file"]
    write_csv(out_csv, rows, fields)
    write_json(out_csv.replace(".csv", "_conversion_report.json"), out)
    log(f"convert-cod: {out}")
    return out
