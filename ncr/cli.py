"""Command-line interface: python -m ncr --config config/study.json <command>.

Steps, in order (docs/REPRODUCE.md):
  convert-oqmd, convert-mp, export-cod, convert-cod, convert-oqmd-icsd, convert-mp-icsd,
  merge-experiments, check, hull, reconstruction, cohorts, power, premise, coverage-known,
  match, estimate, robustness, figures
"""
from __future__ import annotations

import argparse
import glob
import os
import sys
import warnings

warnings.filterwarnings("ignore", message="No Pauling electronegativity")   # noble gases; harmless

from . import __version__                                       # noqa: E402
from . import checks, cohorts, matching, snapshots              # noqa: E402
from . import config as C                                       # noqa: E402
from .util import Manifest, log, read_csv, sha256_file, versions, write_json  # noqa: E402


def _cohort_files(cfg, name):
    files = sorted(glob.glob(os.path.join(cfg["workdir"], "cohorts", "*.csv")))
    if name == "all":
        return files
    p = os.path.join(cfg["workdir"], "cohorts", f"{name}.csv")
    if not os.path.exists(p):
        have = [os.path.splitext(os.path.basename(x))[0] for x in files]
        raise SystemExit(f"no cohort '{name}'; available: {', '.join(have) or 'none (run cohorts first)'}")
    return [p]


def cmd_check(cfg):
    """Validate every input and record its SHA-256 in the manifest (paths relative to the config)."""
    m = Manifest(cfg["workdir"])
    files = {}
    for r in C.releases(cfg):
        p = C.resolve(cfg, r["path"])
        rows = snapshots.read_snapshot(p, r.get("columns"))
        files[r["release_id"]] = {"path": C.relative(cfg, p), "sha256": sha256_file(p), "rows": len(rows),
                                  "hypothetical": sum(1 for x in rows if x["is_hypothetical"])}
        log(f"check: {r['release_id']}: {len(rows)} rows, {files[r['release_id']]['hypothetical']} hypothetical")
    ep = cfg.get("experiments", {}).get("path")
    if ep:
        p = C.resolve(cfg, ep)
        files["experiments"] = {"path": C.relative(cfg, p), "sha256": sha256_file(p),
                                "rows": len(matching.read_experiments(p))}
        log(f"check: experiments: {files['experiments']['rows']} experimental records")
    m.add("inputs", files=files, versions=versions(), pipeline=__version__)
    log(f"check: cutoff {C.cutoff(cfg)}; inputs recorded in {C.relative(cfg, m.path)}")


def cmd_list(cfg):
    for f in _cohort_files(cfg, "all"):
        n = os.path.splitext(os.path.basename(f))[0]
        has_match = os.path.exists(os.path.join(cfg["workdir"], "matches", f"{n}.csv"))
        print(f"{n:30s} rows={len(read_csv(f)):<8d} matched={has_match}")


def build_parser():
    p = argparse.ArgumentParser(prog="python -m ncr", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", default="config/study.json", help="study configuration (JSON)")
    p.add_argument("--version", action="version", version=f"ncr {__version__}")
    sub = p.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("convert-oqmd", help="release snapshot from OQMD MySQL exports (entries and structures TSV)")
    s.add_argument("--entries", required=True); s.add_argument("--structures"); s.add_argument("--out", required=True)
    s.add_argument("--fit", default="standard")
    s = sub.add_parser("convert-mp", help="release snapshot from a Materials Project build (thermo and provenance)")
    s.add_argument("--build-dir", required=True); s.add_argument("--out", required=True)
    s.add_argument("--thermo-type", default="GGA_GGA+U")
    s = sub.add_parser("export-cod", help="export the COD data table from its public server")
    s.add_argument("--out", required=True)
    s = sub.add_parser("convert-cod", help="experimental records from the COD export")
    s.add_argument("--in", dest="inp", required=True); s.add_argument("--out", required=True)
    s = sub.add_parser("convert-oqmd-icsd", help="dated ICSD-derived records from an OQMD export")
    s.add_argument("--in", dest="inp", required=True); s.add_argument("--out", required=True)
    s = sub.add_parser("convert-mp-icsd", help="dated ICSD-derived records from a Materials Project build")
    s.add_argument("--build-dir", required=True); s.add_argument("--out", required=True)
    s = sub.add_parser("merge-experiments", help="merge experimental-record tables into one")
    s.add_argument("--in", dest="inputs", nargs="+", required=True); s.add_argument("--out", required=True)
    sub.add_parser("check", help="validate inputs and record their SHA-256 in the manifest")
    s = sub.add_parser("hull", help="recompute decomposition energies for every release")
    s.add_argument("--workers", type=int, default=1)
    s.add_argument("--all", action="store_true", help="all entries, not only hypothetical ones")
    sub.add_parser("reconstruction", help="compare recomputed and displayed stability labels")
    sub.add_parser("cohorts", help="build the prediction cohorts and the threshold population")
    sub.add_parser("power", help="window counts, minimum detectable effect and check (b)")
    sub.add_parser("premise", help="check (a): cross-database label disagreement")
    sub.add_parser("coverage-known", help="coverage of the outcome record on known compounds")
    s = sub.add_parser("match", help="match predictions against the experimental records")
    s.add_argument("--cohort", default="all", help="cohort name, or 'all'")
    sub.add_parser("estimate", help="all estimates and figures (results/summary.json)")
    sub.add_parser("robustness", help="robustness and sensitivity analyses (results/robustness.json)")
    sub.add_parser("figures", help="Figures 1 to 3 and Supplementary Figures 1 and 2")
    sub.add_parser("list", help="list cohorts and whether they have been matched")
    return p


def main(argv=None):
    a = build_parser().parse_args(argv)
    if a.cmd == "convert-oqmd":
        from .oqmd import convert_oqmd
        return convert_oqmd(a.entries, a.structures, a.out, a.fit)
    if a.cmd == "convert-mp":
        from .mp import convert_mp
        return convert_mp(a.build_dir, a.out, a.thermo_type)
    if a.cmd == "export-cod":
        from .cod import export_cod
        return export_cod(a.out)
    if a.cmd == "convert-cod":
        from .cod import convert_cod
        return convert_cod(a.inp, a.out)
    if a.cmd == "convert-oqmd-icsd":
        from .sources import convert_oqmd_icsd
        return convert_oqmd_icsd(a.inp, a.out)
    if a.cmd == "convert-mp-icsd":
        from .sources import convert_mp_icsd
        return convert_mp_icsd(a.build_dir, a.out)
    if a.cmd == "merge-experiments":
        from .sources import merge_experiments
        return merge_experiments(a.inputs, a.out)
    cfg = C.load(a.config)
    ck = os.path.join(cfg["workdir"], "checks")
    if a.cmd == "check":
        return cmd_check(cfg)
    if a.cmd == "list":
        return cmd_list(cfg)
    if a.cmd == "hull":
        return cohorts.step_hull(cfg, a.workers, a.all)
    if a.cmd == "reconstruction":
        return write_json(os.path.join(ck, "reconstruction.json"), checks.step_reconstruction(cfg))
    if a.cmd == "cohorts":
        return cohorts.step_cohorts(cfg)
    if a.cmd == "power":
        return write_json(os.path.join(ck, "power.json"), checks.step_power(cfg))
    if a.cmd == "premise":
        return write_json(os.path.join(ck, "premise.json"), checks.step_premise(cfg))
    if a.cmd == "coverage-known":
        return write_json(os.path.join(ck, "coverage_known.json"), checks.step_coverage_known(cfg))
    if a.cmd == "match":
        exps = matching.read_experiments(C.resolve(cfg, C.need(cfg, "experiments.path")))
        index = matching.build_exp_index(exps)
        for f in _cohort_files(cfg, a.cohort):
            matching.step_match(cfg, f, exps=exps, index=index)
        return None
    if a.cmd == "estimate":
        from .analysis import step_estimate
        return step_estimate(cfg)
    if a.cmd == "robustness":
        from .robustness import step_robustness
        return step_robustness(cfg)
    if a.cmd == "figures":
        from .analysis import step_figures
        return step_figures(cfg)
    return None


def run():
    try:
        main()
    except (C.ConfigError, FileNotFoundError, RuntimeError, ValueError, KeyError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        if os.environ.get("NCR_DEBUG"):
            raise
        sys.exit(1)
