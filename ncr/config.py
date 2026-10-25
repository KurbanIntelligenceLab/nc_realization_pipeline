"""Configuration: one JSON file describes every input file and every parameter.
Values marked "SET_ME" must be filled before the step that needs them runs."""
from __future__ import annotations

import datetime as dt
import os

from .util import parse_date, read_json

SET_ME = "SET_ME"


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


def need(cfg: dict, dotted: str):
    """Return a config value; raise if missing or still SET_ME."""
    cur = cfg
    for part in dotted.split("."):
        if not isinstance(cur, dict) or part not in cur:
            raise ConfigError(f"config is missing '{dotted}'")
        cur = cur[part]
    if _has_set_me(cur):
        raise ConfigError(f"config value '{dotted}' is still SET_ME; fill it in first")
    return cur


def _has_set_me(v) -> bool:
    if isinstance(v, str):
        return SET_ME in v
    if isinstance(v, dict):
        return any(_has_set_me(x) for x in v.values())
    if isinstance(v, (list, tuple)):
        return any(_has_set_me(x) for x in v)
    return False


def cutoff(cfg: dict) -> dt.date:
    """Data cutoff C: the date of the outcome-record export (Methods). It may not lie in the
    future, because records after the export date cannot have been observed."""
    d, _ = parse_date(need(cfg, "cutoff_date"))
    if d > dt.date.today() and not cfg.get("allow_future_cutoff_for_testing", False):
        raise ConfigError(f"cutoff_date {d} is in the future; use the date of the COD export")
    return d


def releases(cfg: dict, source: str | None = None) -> list[dict]:
    out = []
    if "releases" not in cfg or not cfg["releases"]:
        raise ConfigError("config has no 'releases'")
    for i, r in enumerate(cfg["releases"]):
        for k in ("source", "release_id", "release_date", "path"):
            if k not in r or _has_set_me(r[k]):
                raise ConfigError(f"releases[{i}].{k} is missing or still SET_ME; fill it in first")
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
