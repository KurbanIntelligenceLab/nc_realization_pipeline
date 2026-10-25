#!/usr/bin/env python3
"""Command-line entry point. Run `python3 ncr_cli.py --help` or see RUN_ON_MAC.md.

Typical order (each step writes files under the work directory and can be rerun;
finished outputs are not recomputed unless you delete them):

  check -> hull -> x0 -> cohorts -> [relabel-pop] -> x10 -> premise
  -> for each cohort: dbmatch -> retrieve -> extract -> label -> freeze -> sample
     -> (experts fill sheets) -> adjudicate [--finalize] -> [audit-sample]
  -> estimate
"""
from __future__ import annotations

import argparse
import glob
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import warnings                                        # noqa: E402
warnings.filterwarnings("ignore", message="No Pauling electronegativity")  # noble gases; harmless

from ncr import config as C                     # noqa: E402
from ncr import cohorts, estimate, extract, report, retrieve, sample, snapshots  # noqa: E402
from ncr.util import Manifest, log, read_csv, sha256_file, versions, write_csv, write_json  # noqa: E402


def _cohort_file(cfg, name):
    p = os.path.join(cfg["workdir"], "cohorts", f"{name}.csv")
    if not os.path.exists(p):
        have = [os.path.splitext(os.path.basename(x))[0] for x in glob.glob(os.path.join(cfg["workdir"], "cohorts", "*.csv"))]
        raise SystemExit(f"no cohort '{name}'. Available: {', '.join(sorted(have)) or 'none (run cohorts first)'}")
    return p


def _labels_file(cfg, name, cache="extracted"):
    sub = "labels" if cache == "extracted" else f"labels_{cache}"
    return os.path.join(cfg["workdir"], sub, f"{name}.csv")


def _names(cfg, arg):
    if arg != "all":
        return [arg]
    return sorted(os.path.splitext(os.path.basename(x))[0]
                  for x in glob.glob(os.path.join(cfg["workdir"], "cohorts", "*.csv")))


def cmd_init(a):
    os.makedirs(a.dir, exist_ok=True)
    dst = os.path.join(a.dir, "config.json")
    if os.path.exists(dst):
        raise SystemExit(f"{dst} exists; not overwriting")
    shutil.copy(os.path.join(HERE, "config.example.json"), dst)
    print(f"wrote {dst}; edit it, then run: python3 {os.path.relpath(__file__)} --config {dst} check")


def cmd_check(cfg, a):
    m = Manifest(cfg["workdir"])
    files = {}
    for r in C.releases(cfg):
        p = C.resolve(cfg, r["path"])
        rows = snapshots.read_snapshot(p, r.get("columns"))
        files[r["release_id"]] = {"path": p, "sha256": sha256_file(p), "rows": len(rows),
                                  "hypothetical": sum(1 for x in rows if x["is_hypothetical"])}
        log(f"check: {r['release_id']}: {files[r['release_id']]['rows']} rows, "
            f"{files[r['release_id']]['hypothetical']} hypothetical")
    ep = cfg.get("experiments", {}).get("path")
    if ep and ep != C.SET_ME:
        p = C.resolve(cfg, ep)
        files["experiments"] = {"path": p, "sha256": sha256_file(p), "rows": len(cohorts.read_experiments(p))}
        log(f"check: experiments: {files['experiments']['rows']} experimental records")
    lc = cfg.get("retrieval", {}).get("local_corpus")
    if cfg.get("retrieval", {}).get("backend", "local") == "local" and lc and lc != C.SET_ME:
        p = C.resolve(cfg, lc)
        files["corpus"] = {"path": p, "sha256": sha256_file(p)}
    m.add("inputs", files=files, versions=versions())
    log(f"check: cutoff {C.cutoff(cfg)}; inputs hashed into {m.path}")


def cmd_list(cfg, a):
    for n in _names(cfg, "all"):
        lab = os.path.exists(_labels_file(cfg, n))
        frz = os.path.exists(_labels_file(cfg, n) + ".sha256")
        smp = os.path.exists(os.path.join(cfg["workdir"], "experts", n, "sample_key.csv"))
        fin = os.path.exists(os.path.join(cfg["workdir"], "experts", n, "final_labels.csv"))
        rows = len(read_csv(_cohort_file(cfg, n)))
        print(f"{n:45s} rows={rows:<8d} labeled={lab!s:5s} frozen={frz!s:5s} sampled={smp!s:5s} final={fin!s:5s}")


def cmd_control_cohort(cfg, a):
    rows = read_csv(C.resolve(cfg, C.need(cfg, "control_panel.path")))
    out = []
    for i, r in enumerate(rows):
        for k in ("formula", "panel_kind", "entry_date"):
            if not r.get(k):
                raise SystemExit(f"control panel row {i + 2}: '{k}' is empty")
        f = cohorts.compound_features(r["formula"])
        out.append({"pred_id": f"control::{i}", "cohort": "control_panel", "source": "control",
                    "granularity": r.get("granularity", "composition"), "entry_release": "control",
                    "entry_date": r["entry_date"], "formula": r["formula"], **f,
                    "structure_file": r.get("structure_file", ""), "spacegroup_number": r.get("spacegroup_number", ""),
                    "prototype": "", "panel_kind": r["panel_kind"],
                    "known_first_report_date": r.get("known_first_report_date", "")})
    write_csv(os.path.join(cfg["workdir"], "cohorts", "control_panel.csv"), out)
    log(f"control-cohort: {len(out)} entries")


def cmd_rerun_cohort(cfg, a):
    import numpy as np
    rows = read_csv(_cohort_file(cfg, a.cohort))
    rng = np.random.default_rng(int(cfg["seeds"]["rerun"]))
    idx = sorted(rng.choice(len(rows), size=min(a.n, len(rows)), replace=False).tolist())
    out = os.path.join(cfg["workdir"], "rerun_cohorts", f"{a.cohort}__rerun.csv")
    write_csv(out, [rows[i] for i in idx])
    log(f"rerun-cohort: {len(idx)} entries written to {out}; then run extract/label with --rerun")


def cmd_plan(cfg, a):
    """Count language-model calls before spending anything: one call per retrieved document."""
    from ncr.retrieve import _safe
    from ncr.util import read_json
    total_e = total_d = missing = 0
    for name in _names(cfg, a.cohort):
        rows = read_csv(_cohort_file(cfg, name))
        docs = 0
        for r in rows:
            p = os.path.join(cfg["workdir"], "retrieved", _safe(r["pred_id"]) + ".json")
            if os.path.exists(p):
                docs += len(read_json(p)["docs"])
            else:
                missing += 1
        total_e += len(rows)
        total_d += docs
        print(f"{name:45s} entries={len(rows):<8d} documents={docs}")
    print(f"total: {total_e} entries, {total_d} documents = {total_d} model calls"
          + (f"; {missing} entries not yet retrieved" if missing else ""))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", default="config.json")
    sub = p.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("init", help="write a config template"); s.add_argument("--dir", default=".")
    sub.add_parser("check", help="validate inputs and record their SHA-256")
    sub.add_parser("list", help="show cohorts and progress")
    s = sub.add_parser("convert-entries", help="snapshot CSV from a pymatgen entries JSON")
    s.add_argument("--in", dest="inp", required=True); s.add_argument("--out", required=True)
    s.add_argument("--hypothetical-ids"); s.add_argument("--all-hypothetical", action="store_true")
    s.add_argument("--cif-dir"); s.add_argument("--id-field", default="entry_id")
    s = sub.add_parser("convert-oqmd", help="snapshot CSV from OQMD MySQL exports (entries and structures TSV)")
    s.add_argument("--entries", required=True); s.add_argument("--structures"); s.add_argument("--out", required=True)
    s.add_argument("--fit", default="standard")
    s = sub.add_parser("export-cod", help="export the COD data table from its public server (no CIFs)")
    s.add_argument("--out", required=True)
    s = sub.add_parser("convert-cod", help="experimental-reports table from the COD export")
    s.add_argument("--in", dest="inp", required=True); s.add_argument("--out", required=True)
    s = sub.add_parser("convert-oqmd-icsd", help="dated ICSD-derived records from an OQMD export (entry_id, composition_id, year)")
    s.add_argument("--in", dest="inp", required=True); s.add_argument("--out", required=True)
    s = sub.add_parser("convert-mp-icsd", help="dated ICSD-derived records from a Materials Project build (provenance + thermo)")
    s.add_argument("--build-dir", required=True); s.add_argument("--out", required=True)
    s = sub.add_parser("merge-experiments", help="merge experimental-report tables into one")
    s.add_argument("--in", dest="inputs", nargs="+", required=True); s.add_argument("--out", required=True)
    s = sub.add_parser("convert-mp", help="snapshot CSV from a Materials Project build (thermo + provenance)")
    s.add_argument("--build-dir", required=True); s.add_argument("--out", required=True)
    s.add_argument("--thermo-type", default="GGA_GGA+U")
    s = sub.add_parser("hull", help="recompute decomposition energies for every release")
    s.add_argument("--workers", type=int, default=1)
    s.add_argument("--all", action="store_true", help="all entries, not only hypothetical (needed for 'premise')")
    sub.add_parser("x0", help="snapshot reconstruction counts")
    sub.add_parser("cohorts", help="build prediction cohorts and threshold populations")
    sub.add_parser("relabel-pop", help="build relabeling populations (correction-driven flips)")
    sub.add_parser("x10", help="window census and minimum detectable effect")
    sub.add_parser("premise", help="E6a cross-database label disagreement")
    sub.add_parser("control-cohort", help="build the X1 control-panel cohort")
    sub.add_parser("coverage-known", help="design-stage check: known (ICSD-origin) compositions found in the reports table")
    for name, hlp in (("dbmatch", "match a cohort against the structure-database export"),
                      ("retrieve", "literature retrieval"), ("extract", "language-model extraction"),
                      ("label", "entry-level model labels"), ("freeze", "hash and lock the labels"),
                      ("sample", "draw the expert sample"), ("adjudicate", "merge expert sheets"),
                      ("audit-sample", "X3 residual-miss sample")):
        s = sub.add_parser(name, help=hlp)
        s.add_argument("--cohort", required=True, help="cohort name, or 'all'")
        if name in ("extract",):
            s.add_argument("--provider", choices=["anthropic", "mock"], default="anthropic")
        if name in ("extract", "label", "retrieve"):
            s.add_argument("--rerun", action="store_true", help="use the X5 rerun subset and cache")
        if name == "sample":
            s.add_argument("--n", type=int)
        if name == "adjudicate":
            s.add_argument("--finalize", action="store_true")
    s = sub.add_parser("plan", help="count model calls before running extract")
    s.add_argument("--cohort", default="all")
    s = sub.add_parser("rerun-cohort", help="X5: random subset for a second extraction run")
    s.add_argument("--cohort", required=True); s.add_argument("--n", type=int, default=300)
    sub.add_parser("estimate-census", help="annotation-free estimates from dbmatch outputs (deterministic design)")
    sub.add_parser("robustness", help="all robustness and sensitivity analyses (Supplementary Note 5)")
    s = sub.add_parser("estimate", help="all certified estimates, figures and results_for_claude.json")
    s.add_argument("--synthetic", action="store_true", help=argparse.SUPPRESS)
    sub.add_parser("version")
    a = p.parse_args(argv)

    if a.cmd == "init":
        return cmd_init(a)
    if a.cmd == "version":
        print(versions())
        return
    if a.cmd == "convert-entries":
        n = snapshots.convert_entries(a.inp, a.out, a.hypothetical_ids, a.all_hypothetical, a.cif_dir, a.id_field)
        log(f"convert-entries: {n} entries written to {a.out}")
        return
    if a.cmd == "export-cod":
        from ncr.cod import export_cod
        export_cod(a.out)
        return
    if a.cmd == "convert-cod":
        from ncr.cod import convert_cod
        convert_cod(a.inp, a.out)
        return
    if a.cmd == "convert-oqmd-icsd":
        from ncr.sources import convert_oqmd_icsd
        convert_oqmd_icsd(a.inp, a.out); return
    if a.cmd == "convert-mp-icsd":
        from ncr.sources import convert_mp_icsd
        convert_mp_icsd(a.build_dir, a.out); return
    if a.cmd == "merge-experiments":
        from ncr.sources import merge_experiments
        merge_experiments(a.inputs, a.out); return
    if a.cmd == "convert-mp":
        from ncr.mp import convert_mp
        convert_mp(a.build_dir, a.out, a.thermo_type)
        return
    if a.cmd == "convert-oqmd":
        from ncr.oqmd import convert_oqmd
        convert_oqmd(a.entries, a.structures, a.out, a.fit)
        return
    cfg = C.load(a.config)
    m = Manifest(cfg["workdir"])
    if a.cmd == "check":
        return cmd_check(cfg, a)
    if a.cmd == "list":
        return cmd_list(cfg, a)
    if a.cmd == "hull":
        return cohorts.step_hull(cfg, a.workers, a.all)
    if a.cmd == "x0":
        write_json(os.path.join(cfg["workdir"], "checks", "x0.json"), cohorts.step_x0(cfg)); return
    if a.cmd == "cohorts":
        return cohorts.step_cohorts(cfg)
    if a.cmd == "relabel-pop":
        return cohorts.step_relabel(cfg)
    if a.cmd == "x10":
        write_json(os.path.join(cfg["workdir"], "checks", "x10.json"), estimate.step_x10(cfg)); return
    if a.cmd == "premise":
        write_json(os.path.join(cfg["workdir"], "checks", "premise.json"), estimate.step_premise(cfg)); return
    if a.cmd == "coverage-known":
        write_json(os.path.join(cfg["workdir"], "checks", "coverage_known.json"), estimate.step_coverage_known(cfg)); return
    if a.cmd == "control-cohort":
        return cmd_control_cohort(cfg, a)
    if a.cmd == "rerun-cohort":
        return cmd_rerun_cohort(cfg, a)
    if a.cmd == "plan":
        return cmd_plan(cfg, a)
    if a.cmd == "robustness":
        from ncr.robust import step_robustness
        step_robustness(cfg); return
    if a.cmd == "estimate-census":
        from ncr.census import step_estimate_census
        step_estimate_census(cfg); return
    if a.cmd == "estimate":
        report.step_estimate(cfg, synthetic=a.synthetic); return

    cache = "extracted_rerun" if getattr(a, "rerun", False) else "extracted"
    for name in _names(cfg, a.cohort):
        if getattr(a, "rerun", False):
            cf = os.path.join(cfg["workdir"], "rerun_cohorts", f"{name}__rerun.csv")
        else:
            cf = _cohort_file(cfg, name)
        if a.cmd == "dbmatch":
            cohorts.step_dbmatch(cfg, cf)
        elif a.cmd == "retrieve":
            retrieve.step_retrieve(cfg, cf, lambda x: C.resolve(cfg, x), str(C.cutoff(cfg)))
        elif a.cmd == "extract":
            extract.step_extract(cfg, cf, a.provider, cache)
        elif a.cmd == "label":
            if getattr(a, "rerun", False):
                extract.step_label(cfg, cf, cache,
                                   out_path=os.path.join(cfg["workdir"], "labels_extracted_rerun", f"{name}__rerun.csv"),
                                   db_path=os.path.join(cfg["workdir"], "dbmatch", f"{name}.csv"))
            else:
                extract.step_label(cfg, cf, cache)
        elif a.cmd == "freeze":
            extract.step_freeze(cfg, _labels_file(cfg, name), m)
        elif a.cmd == "sample":
            sample.step_sample(cfg, cf, _labels_file(cfg, name), m, a.n)
        elif a.cmd == "adjudicate":
            sample.step_adjudicate(cfg, cf, _labels_file(cfg, name), a.finalize)
        elif a.cmd == "audit-sample":
            sample.step_audit_sample(cfg, cf, _labels_file(cfg, name))


if __name__ == "__main__":
    try:
        main()
    except (C.ConfigError, FileNotFoundError, RuntimeError, ValueError, KeyError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        if os.environ.get("NCR_DEBUG"):
            raise
        sys.exit(1)
