"""Study configuration: one JSON file lists every input file and every parameter.
Relative paths are resolved against the directory that contains the file."""
from __future__ import annotations

import datetime as dt
import os

from .util import parse_date, read_json


class ConfigError(ValueError):
    pass


def load(path: str) -> dict:
    cfg = read_json(path)
    cfg["_dir"] = os.path.dirname(os.path.abspath(path))
    cfg.setdefault("workdir", "work")
    cfg["workdir"] = resolve(cfg, cfg["workdir"])
    os.makedirs(cfg["workdir"], exist_ok=True)
    return cfg


def resolve(cfg: dict, p: str) -> str:
    return p if os.path.isabs(p) else os.path.join(cfg["_dir"], p)


def relative(cfg: dict, p: str) -> str:
    """Path relative to the configuration directory, for portable records."""
    return os.path.relpath(os.path.abspath(p), cfg["_dir"])


def need(cfg: dict, dotted: str):
    """Return a configuration value; raise ConfigError if it is missing."""
    cur = cfg
    for part in dotted.split("."):
        if not isinstance(cur, dict) or part not in cur:
            raise ConfigError(f"configuration is missing '{dotted}'")
        cur = cur[part]
    return cur


def cutoff(cfg: dict) -> dt.date:
    """Data cutoff C: the date of the outcome-record export. It may not lie in the future,
    because records after the export date cannot have been observed."""
    d, _ = parse_date(need(cfg, "cutoff_date"))
    if d > dt.date.today() and not cfg.get("allow_future_cutoff", False):
        raise ConfigError(f"cutoff_date {d} is in the future; use the date of the COD export")
    return d


def releases(cfg: dict, source: str | None = None) -> list[dict]:
    if not cfg.get("releases"):
        raise ConfigError("configuration has no 'releases'")
    out = []
    for i, r in enumerate(cfg["releases"]):
        for k in ("source", "release_id", "release_date", "path"):
            if k not in r:
                raise ConfigError(f"releases[{i}] is missing '{k}'")
        if source is None or r["source"] == source:
            r = dict(r)
            r["date"], _ = parse_date(r["release_date"])
            out.append(r)
    return sorted(out, key=lambda r: r["date"])


def sources(cfg: dict) -> list[str]:
    seen = []
    for r in releases(cfg):
        if r["source"] not in seen:
            seen.append(r["source"])
    return seen
