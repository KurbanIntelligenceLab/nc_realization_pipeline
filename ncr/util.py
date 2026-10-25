"""Shared helpers: dates, hashing, CSV I/O, logging and the run manifest."""
from __future__ import annotations

import csv
import datetime as dt
import hashlib
import json
import os
import sys
import time

YEAR_ONLY_MONTH, YEAR_ONLY_DAY = 7, 1      # Table S1: year-only dates set to 1 July
MONTH_ONLY_DAY = 15                        # implementation choice (README, "Choices to confirm")


def parse_date(s):
    """Parse YYYY-MM-DD, YYYY-MM or YYYY. Returns (datetime.date, precision) or
    (None, None) for empty values or 'none'."""
    if s is None:
        return None, None
    s = str(s).strip()
    if s == "" or s.lower() in {"none", "nan", "null", "na", "n/a"}:
        return None, None
    parts = s.split("T")[0].split("-")
    try:
        if len(parts) == 1:
            return dt.date(int(parts[0]), YEAR_ONLY_MONTH, YEAR_ONLY_DAY), "year"
        if len(parts) == 2:
            return dt.date(int(parts[0]), int(parts[1]), MONTH_ONLY_DAY), "month"
        return dt.date(int(parts[0]), int(parts[1]), int(parts[2])), "day"
    except ValueError as exc:
        raise ValueError(f"unreadable date {s!r}; use YYYY-MM-DD, YYYY-MM or YYYY") from exc


def decimal_year(d: dt.date | None) -> float:
    if d is None:
        return float("nan")
    start = dt.date(d.year, 1, 1)
    days = (dt.date(d.year + 1, 1, 1) - start).days
    return d.year + (d - start).days / days


def date_str(d: dt.date | None) -> str:
    return "" if d is None else d.isoformat()


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_csv(path: str) -> list[dict]:
    if not os.path.exists(path):
        raise FileNotFoundError(f"missing input file: {path}")
    with open(path, newline="", encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh))


def write_csv(path: str, rows: list[dict], fields: list[str] | None = None) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    if fields is None:
        fields = []
        for r in rows:
            for k in r:
                if k not in fields:
                    fields.append(k)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: ("" if r.get(k) is None else r.get(k)) for k in fields})


def read_json(path: str):
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def write_json(path: str, obj) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(obj, fh, indent=2, default=_json_default)


def _json_default(o):
    try:
        import numpy as np
        if isinstance(o, np.bool_):
            return bool(o)
        if isinstance(o, np.integer):
            return int(o)
        if isinstance(o, np.floating):
            return float(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
    except ImportError:  # pragma: no cover
        pass
    if isinstance(o, (dt.date, dt.datetime)):
        return o.isoformat()
    if isinstance(o, (set, tuple)):
        return list(o)
    raise TypeError(f"not JSON serializable: {type(o)}")


def truthy(v) -> bool:
    return str(v).strip().lower() in {"1", "true", "yes", "y", "t"}


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", file=sys.stderr, flush=True)


class Manifest:
    """Append-only record of inputs, hashes, seeds and versions (work/manifest.json)."""

    def __init__(self, workdir: str):
        self.path = os.path.join(workdir, "manifest.json")
        self.data = read_json(self.path) if os.path.exists(self.path) else {"events": []}

    def add(self, kind: str, **kw) -> None:
        kw.update({"kind": kind, "time": dt.datetime.now().isoformat(timespec="seconds")})
        self.data["events"].append(kw)
        write_json(self.path, self.data)

    def last(self, kind: str, **match):
        for ev in reversed(self.data["events"]):
            if ev.get("kind") == kind and all(ev.get(k) == v for k, v in match.items()):
                return ev
        return None


def versions() -> dict:
    out = {"python": sys.version.split()[0]}
    for mod in ("numpy", "pymatgen", "matplotlib", "anthropic"):
        try:
            from importlib.metadata import version
            out[mod] = version(mod)
        except Exception:  # noqa: BLE001
            out[mod] = None
    return out
