"""Run every estimate that has its inputs, write the full results, the figure JSON,
the figures, formatted table strings, and results_for_claude.json (aggregates only)."""
from __future__ import annotations

import glob
import os
import sys

import numpy as np

from . import config as C
from .estimate import (Bundle, _clip, _y, base_rates, control_metrics, relabel, rerun_compare,
                       residual_bound, threshold)
from .util import Manifest, log, read_csv, read_json, versions, write_json

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)
import nc_realization_metrics as NCM                 # noqa: E402


def _fmt(d):
    try:
        return NCM.format_share(d["est"], min(d["lo"], d["est"]), max(d["hi"], d["est"]))
    except Exception:  # noqa: BLE001
        return None


def step_estimate(cfg: dict, synthetic: bool = False) -> dict:
    wd = cfg["workdir"]
    res = {"base": {}, "threshold": {}, "relabel": {}, "status": {}}
    for cf in sorted(glob.glob(os.path.join(wd, "cohorts", "*.csv"))):
        stem = os.path.splitext(os.path.basename(cf))[0]
        lf = os.path.join(wd, "labels", stem + ".csv")
        fin = os.path.join(wd, "experts", stem, "final_labels.csv")
        if stem == "control_panel":
            if os.path.exists(lf):
                res["control_panel"] = control_metrics(cfg, cf, lf)
            continue
        if not os.path.exists(fin):
            res["status"][stem] = "no final expert labels yet"
            continue
        b = Bundle(cfg, cf, lf)
        if stem.startswith("threshold__"):
            res["threshold"][stem] = threshold(cfg, b)
        elif stem.startswith("relabel__"):
            res["relabel"][stem] = relabel(cfg, b)
        else:
            full = cf.replace("__phase1.csv", ".csv")
            n_pop = len(read_csv(full)) if stem.endswith("__phase1") else None
            dbp = os.path.join(wd, "dbmatch", stem + ".csv")
            if os.path.exists(dbp):
                gran = b.R[0].get("granularity", "structure")
                col = "db_first_T2" if gran == "structure" else "db_first_T3"
                db = {x["pred_id"]: _y(x[col]) for x in read_csv(dbp)}
                cfg.setdefault("_db_only", {})[stem] = np.array([db.get(x["pred_id"], np.nan) for x in b.labs])
            src = b.R[0]["source"]
            r = base_rates(cfg, b, src, n_pop)
            main_h = str(cfg.get("ml_horizons_years", [2])[0] if src in cfg.get("ml_sources", [])
                         else int(cfg.get("threshold", {}).get("horizon_years", 5)))
            r["residual_miss"] = residual_bound(cfg, b, r, main_h)
            res["base"][stem] = r
        res["status"][stem] = "estimated"
    # Rerun agreement (X5)
    res["rerun"] = {}
    for lf in glob.glob(os.path.join(wd, "labels_extracted_rerun", "*.csv")):
        stem = os.path.splitext(os.path.basename(lf))[0]
        main = os.path.join(wd, "labels", stem.replace("__rerun", "") + ".csv")
        if os.path.exists(main):
            res["rerun"][stem] = rerun_compare(main, lf)
    for name in ("x0", "x10", "premise"):
        p = os.path.join(wd, "checks", f"{name}.json")
        if os.path.exists(p):
            res[name] = read_json(p)
    # Figures
    figs = _figure_json(cfg, res)
    if synthetic:
        figs["SYNTHETIC_demo"] = True
    res["figures"] = figs
    if figs.get("fig2a") and figs.get("fig3a"):
        try:
            import nc_realization_figures as NCF
            res["figure_checks"] = NCF.render_all(figs, os.path.join(wd, "figures"), synthetic=synthetic)
        except Exception as exc:  # noqa: BLE001
            res["figure_checks"] = {"error": str(exc)}
    res["tables"] = _tables(res)
    m = Manifest(wd)
    res["provenance"] = {"versions": versions(), "freezes": [e for e in m.data["events"] if e["kind"] == "freeze"],
                         "samples": [{k: e[k] for k in ("labels_file", "labels_sha256", "seed", "n")}
                                     for e in m.data["events"] if e["kind"] == "expert_sample"],
                         "inputs": [e for e in m.data["events"] if e["kind"] == "inputs"]}
    cl = os.path.join(wd, "cost_log.csv")
    res["cost_log"] = read_csv(cl) if os.path.exists(cl) else []
    write_json(os.path.join(wd, "results", "results_full.json"), res)
    write_json(os.path.join(wd, "results_for_claude.json"), _for_claude(cfg, res))
    log(f"estimate: wrote {os.path.join(wd, 'results_for_claude.json')}")
    return res


def _figure_json(cfg, res) -> dict:
    a, b = [], []
    for stem, r in res["base"].items():
        hs = [h for h, v in r["horizons"].items() if "realized_share" in v]
        if not hs:
            continue
        hmax = max(hs, key=float)
        v = r["horizons"][hmax]
        a.append({"source": r["source"], "horizon": float(hmax),
                  "already_known": _clip(v["pre_existence_within_R"]), "realized": _clip(v["realized_share"]),
                  "not_realized": _clip(v["not_realized_share"])})
        for h in hs:
            rr = r["horizons"][h].get("realized_among_not_known", {})
            if "est" in rr:
                b.append({"source": r["source"], "t": float(h), **_clip(rr)})
    out = {"fig2a": a, "fig2b": b}
    src = cfg.get("threshold", {}).get("figure_source")
    for stem, t in res["threshold"].items():
        if src is None or stem.endswith(_slug(src)):
            out["fig3a"] = t["fig3a"]
            break
    return out


def _slug(s):
    return "".join(ch if ch.isalnum() else "_" for ch in s).strip("_").lower()


def _tables(res) -> dict:
    t = {}
    for stem, r in res["base"].items():
        row = {"pre_existence_all": _fmt(r["pre_existence_all"])}
        for h, v in r["horizons"].items():
            if "realized_share" in v:
                row[f"realized_{h}y"] = _fmt(v["realized_share"])
                if "est" in v.get("realized_among_not_known", {}):
                    row[f"realized_among_not_known_{h}y"] = _fmt(v["realized_among_not_known"])
        t[stem] = row
    return t


def _for_claude(cfg, res) -> dict:
    """Aggregates only: no formulas, identifiers, documents or database content."""
    def strip(o):
        if isinstance(o, dict):
            return {k: strip(v) for k, v in o.items() if k not in {"blocks", "bins_detail"}}
        if isinstance(o, list):
            return [strip(v) for v in o]
        return o
    keep = {k: res.get(k) for k in ("base", "threshold", "relabel", "control_panel", "rerun", "x0", "x10",
                                    "premise", "tables", "figure_checks", "status")}
    keep["figures"] = res.get("figures")
    keep["provenance"] = {"versions": res["provenance"]["versions"],
                          "frozen_label_hashes": [e.get("sha256") for e in res["provenance"]["freezes"]],
                          "samples": res["provenance"]["samples"]}
    keep["cost_totals"] = {"calls": sum(int(c.get("calls", 0) or 0) for c in res.get("cost_log", [])),
                           "tokens_in": sum(int(c.get("tokens_in", 0) or 0) for c in res.get("cost_log", [])),
                           "tokens_out": sum(int(c.get("tokens_out", 0) or 0) for c in res.get("cost_log", []))}
    keep["config_summary"] = {k: cfg.get(k) for k in ("cutoff_date", "horizons_years", "ml_horizons_years",
                                                      "ml_sources", "sampling", "threshold", "two_phase")}
    keep["config_summary"]["llm_model"] = cfg.get("llm", {}).get("model")
    return strip(keep)
