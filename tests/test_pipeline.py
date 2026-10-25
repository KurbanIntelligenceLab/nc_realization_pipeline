#!/usr/bin/env python3
"""End-to-end test on a synthetic world (no real data, mock language model,
simulated experts). Run: python3 tests/test_pipeline.py
Checks that every step runs, that the safeguards fire, and that certified
estimates recover the synthetic truth within their uncertainty."""
from __future__ import annotations

import datetime as dt
import glob
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "tests"))
import make_synthetic                                   # noqa: E402
from ncr.util import decimal_year, parse_date, read_csv, read_json, write_csv  # noqa: E402

CLI = [sys.executable, os.path.join(HERE, "ncr_cli.py")]
RESULTS = []


def check(name, cond, detail=""):
    RESULTS.append(bool(cond))
    print(f"[{'PASS' if cond else 'FAIL'}] {name}" + (f": {detail}" if detail else ""))


def run(d, *args, ok=True):
    p = subprocess.run(CLI + ["--config", os.path.join(d, "config.json"), *args], capture_output=True, text=True)
    if ok and p.returncode != 0:
        print(p.stdout, p.stderr)
        raise SystemExit(f"step failed: {' '.join(args)}")
    return p


def simulate_experts(d, truth, rng):
    tmap = {t["formula"]: t for t in truth}
    for sd in glob.glob(os.path.join(d, "work", "experts", "*")):
        for who, miss, jitter in (("A", 0.05, 0.03), ("B", 0.05, 0.03)):
            p = os.path.join(sd, f"expert_{who}.csv")
            rows = read_csv(p)
            for r in rows:
                t = tmap[r["formula"]]
                gran = r["granularity"]
                rep = t["first_report"]
                if rep is not None and rng.random() < miss:
                    rep = None
                if rep is not None and rng.random() < jitter:
                    rep = rep + rng.choice([-1, 1]) * rng.uniform(1.5, 3)
                r["tier"] = ("T2" if gran == "structure" else "T3") if rep is not None else "none"
                r["first_report_date"] = make_synthetic._date(rep) if rep is not None else ""
                r["first_mention_date"] = make_synthetic._date(t["first_mention"]) if t["first_mention"] is not None else ""
                r["minutes_spent"] = "10"
            write_csv(p, rows)


def adjudicate_third(d, truth):
    tmap = {t["formula"]: t for t in truth}
    for p in glob.glob(os.path.join(d, "work", "experts", "*", "adjudication.csv")):
        rows = read_csv(p)
        for r in rows:
            t = tmap[r["formula"]]
            r["final_has_report"] = "yes" if t["first_report"] is not None else "no"
            r["final_first_report_date"] = make_synthetic._date(t["first_report"]) if t["first_report"] is not None else ""
            r["final_first_mention_date"] = make_synthetic._date(t["first_mention"]) if t["first_mention"] is not None else ""
        write_csv(p, rows)


def true_share(d, cohort, truth, fn):
    tmap = {t["formula"]: t for t in truth}
    rows = read_csv(os.path.join(d, "work", "cohorts", f"{cohort}.csv"))
    vals = [fn(tmap[r["formula"]], decimal_year(parse_date(r["entry_date"])[0]), r) for r in rows]
    vals = [v for v in vals if v is not None]
    return float(np.mean(vals)), len(vals)


def main():
    d = tempfile.mkdtemp(prefix="ncr_test_")
    try:
        cfg = make_synthetic.make(d)
        truth = read_json(os.path.join(d, "truth.json"))["truth"]
        rng = np.random.default_rng(99)
        run(d, "check")
        run(d, "hull", "--all")
        run(d, "x0")
        run(d, "cohorts")
        run(d, "x10")
        run(d, "premise")
        x10 = read_json(os.path.join(d, "work", "checks", "x10.json"))
        cs = x10["SrcB"].get("block_scheme_candidates_25meV", {}).get("schemes", {})
        check("x10 compares block schemes (blocked SMD, usable counts, MDE)",
              set(cs) == {"class", "class_en", "class_en_nel"} and all("max_abs_blocked_smd" in v for v in cs.values())
              and set(x10["SrcB"].get("block_scheme_balance_by_window", {})) == {"10", "25", "50"},
              f"{ {k: (round(v['max_abs_blocked_smd'], 3) if v['max_abs_blocked_smd'] is not None else None, v['N_A_usable']) for k, v in cs.items()} }")
        check("x10 balance (check b) computed before outcomes", "check_b_passes" in x10["SrcB"]["25"]["balance"],
              f"max |SMD| {x10['SrcB']['25']['balance']['max_abs_smd_overall']}")
        check("x10 window census", x10["SrcB"]["25"]["N_A"] > 0 and x10["SrcB"]["25"]["N_B"] > 0,
              f"SrcB 25 meV: N_A={x10['SrcB']['25']['N_A']}, N_B={x10['SrcB']['25']['N_B']}")
        prem = read_json(os.path.join(d, "work", "checks", "premise.json"))
        check("premise check runs", "srcb vs srca_r1" in prem and prem["srcb vs srca_r1"]["matched"] > 0)
        names = [os.path.splitext(os.path.basename(x))[0] for x in glob.glob(os.path.join(d, "work", "cohorts", "*.csv"))]
        check("cohorts built", {"srca", "srcb", "threshold__srca", "threshold__srcb"} <= set(names), str(sorted(names)))
        for step in ("dbmatch", "retrieve"):
            run(d, step, "--cohort", "all")
        run(d, "extract", "--cohort", "all", "--provider", "mock")
        run(d, "label", "--cohort", "all")
        # Safeguard: sampling before freezing must fail
        p = run(d, "sample", "--cohort", "srcb", ok=False)
        check("sampling refused before freezing", p.returncode != 0 and "not frozen" in (p.stderr + p.stdout))
        run(d, "freeze", "--cohort", "all")
        # Safeguard: relabeling after freezing must fail
        p = run(d, "label", "--cohort", "srcb", ok=False)
        check("labels locked after freezing", p.returncode != 0)
        run(d, "sample", "--cohort", "all")
        key = read_csv(os.path.join(d, "work", "experts", "srcb", "sample_key.csv"))
        sheet = read_csv(os.path.join(d, "work", "experts", "srcb", "expert_A.csv"))
        blind = set(sheet[0]) & {"pred_id", "model_label", "side", "x_mev", "source", "dhd_entry", "stratum"}
        check("expert sheets are blinded", not blind, f"leaked columns: {blind}" if blind else "no label, side or source")
        simulate_experts(d, truth, rng)
        run(d, "adjudicate", "--cohort", "all")
        adjudicate_third(d, truth)
        run(d, "adjudicate", "--cohort", "all", "--finalize")
        kap = read_json(os.path.join(d, "work", "experts", "srcb", "kappa.json"))
        check("kappa computed", 0.5 < kap["has_report"] <= 1.0, f"{kap}")
        run(d, "audit-sample", "--cohort", "srcb")
        sh = read_csv(os.path.join(d, "work", "experts", "srcb", "residual_audit_sheet.csv"))
        for r in sh:
            r["found_missed_report"] = "no"
        write_csv(os.path.join(d, "work", "experts", "srcb", "residual_audit_sheet.csv"), sh)
        run(d, "control-cohort")
        for step in ("retrieve", "extract", "label"):
            run(d, step, "--cohort", "control_panel", *(["--provider", "mock"] if step == "extract" else []))
        run(d, "rerun-cohort", "--cohort", "srcb", "--n", "60")
        run(d, "extract", "--cohort", "srcb", "--provider", "mock", "--rerun")
        run(d, "label", "--cohort", "srcb", "--rerun")
        run(d, "estimate", "--synthetic")
        res = read_json(os.path.join(d, "work", "results", "results_full.json"))
        rc = read_json(os.path.join(d, "work", "results_for_claude.json"))
        check("results_for_claude written", "base" in rc and "threshold" in rc)
        txt = json.dumps(rc)
        leaked = [t["formula"] for t in truth[:200] if f'"{t["formula"]}"' in txt]
        check("results_for_claude has no formulas", not leaked, f"{leaked[:3]}")

        # Truth recovery: pre-existence (whole cohort) and 5-year realization, per base cohort
        for coh in ("srcb", "srca"):
            r = res["base"][coh]
            pre = r["pre_existence_all"]
            tv, n = true_share(d, coh, truth, lambda t, e, row: float(t["first_report"] is not None and t["first_report"] <= e))
            ok = abs(pre["est"] - tv) <= 4 * pre["se"] + 1e-9
            check(f"{coh}: certified pre-existence recovers truth", ok,
                  f"est {pre['est']:.3f} [{pre['lo']:.3f}, {pre['hi']:.3f}] vs truth {tv:.3f} (N={n})")
            h5 = r["horizons"].get("5", {})
            if "realized_share" in h5:
                cut = decimal_year(dt.date(2026, 12, 31)) - 1
                tv5, n5 = true_share(d, coh, truth, lambda t, e, row: (float(t["first_report"] is not None and e < t["first_report"] <= e + 5)
                                                                     if e + 5 <= cut else None))
                rs = h5["realized_share"]
                check(f"{coh}: certified 5-year realization recovers truth",
                      abs(rs["est"] - tv5) <= 4 * rs["se"] + 1e-9,
                      f"est {rs['est']:.3f} [{rs['lo']:.3f}, {rs['hi']:.3f}] vs truth {tv5:.3f}; model-only {h5['naive_model_only']:.3f}")
        # Threshold effect (SrcB, 25 meV/atom) against the finite-population truth
        th = res["threshold"]["threshold__srcb"]["windows"]["25"]
        rows = read_csv(os.path.join(d, "work", "cohorts", "threshold__srcb.csv"))
        tmap = {t["formula"]: t for t in truth}
        ya, yb = [], []
        for r in rows:
            if abs(float(r["x_mev"])) > 25:
                continue
            e = decimal_year(parse_date(r["entry_date"])[0])
            t = tmap[r["formula"]]["first_report"]
            y = float(t is not None and e < t <= e + 5)
            (ya if r["side"] == "A" else yb).append(y)
        tdiff = np.mean(ya) - np.mean(yb)
        eff = th.get("effect", {})
        check("threshold effect estimated", "estimate" in eff, th.get("note", ""))
        if "estimate" in eff:
            check("threshold effect recovers the finite-population difference",
                  abs(eff["estimate"] - tdiff) <= 4 * eff["se"],
                  f"est {eff['estimate']:.3f} [{eff['ci_low']:.3f}, {eff['ci_high']:.3f}] vs truth {tdiff:.3f}; MDE {th['mde']:.3f}")
        check("balance, density and placebo reported",
              all(k in res["threshold"]["threshold__srcb"] for k in ("balance_primary_window", "density_test", "placebo")))
        check("control panel metrics", res.get("control_panel", {}).get("sensitivity") is not None,
              f"{res.get('control_panel')}")
        check("rerun agreement (X5)", any(v["label_agreement"] == 1.0 for v in res["rerun"].values()),
              f"{res['rerun']}")
        check("residual-miss bound (X3)", "adjusted_upper_bound" in (res["base"]["srcb"].get("residual_miss") or {}),
              f"{res['base']['srcb'].get('residual_miss')}")
        fc = res.get("figure_checks", {})
        check("figures rendered with checks", fc and all(v is True for v in fc.values()), f"{fc}")
        # Relabeling estimator on a synthetic relabeling population (estimation path only)
        relabel_estimation_test(d, truth)
        tools_test(d)
        oqmd_converter_test(d)
        mp_converter_test(d)
        polymorph_test()
        density_test_large_n()
        cod_converter_test(d)
        census_test(d, truth)
        two_phase_ratio_coverage()
    finally:
        if os.environ.get("NCR_KEEP_TEST_DIR"):
            print(f"test directory kept: {d}")
        else:
            shutil.rmtree(d, ignore_errors=True)
    print(f"\n{sum(RESULTS)}/{len(RESULTS)} checks passed")
    if not all(RESULTS):
        raise SystemExit(1)
    print("ALL PIPELINE CHECKS PASSED")


def relabel_estimation_test(d, truth):
    """Estimation path of the relabeling design. The scheme-application step needs real
    Materials Project entries with calculation parameters and is not exercised here."""
    from ncr import estimate as E
    import ncr.config as C
    src = os.path.join(d, "work", "cohorts", "srca.csv")
    rows = read_csv(src)
    rng = np.random.default_rng(5)
    for r in rows:
        r["group"] = "F" if rng.random() < 0.3 else "M"
        r["release_r_date"] = "2014-06-01"
        r["cem_key"] = f"{r['n_elements']}|{r['chemistry_class']}"
    # reuse the srca labels and expert sample under a relabel-type bundle
    cfg = C.load(os.path.join(d, "config.json"))
    b = E.Bundle(cfg, src, os.path.join(d, "work", "labels", "srca.csv"))
    b.R = [{**x, **{k: r[k] for k in ("group", "release_r_date", "cem_key")}}
           for x, r in zip(b.R, [next(rr for rr in rows if rr["pred_id"] == l["pred_id"]) for l in b.labs])]
    out = E.relabel(cfg, b)
    check("relabeling estimator runs (treated weights, CEM blocks)",
          "estimate" in out.get("effect", {}) and "estimate" in out.get("pre_period_placebo", {}),
          f"counts {out['counts']}; effect {out.get('effect', {}).get('estimate')}")


def tools_test(d):
    """convert-entries, parallel hull and the relabel-pop scheme step (smoke test)."""
    from pymatgen.entries.computed_entries import ComputedEntry
    from ncr.snapshots import compute_dhd, read_snapshot
    ents = [ComputedEntry("Li", -1.9, entry_id="e1").as_dict(), ComputedEntry("O2", -9.8, entry_id="e2").as_dict(),
            ComputedEntry("Li2O", -14.3, entry_id="e3").as_dict(), ComputedEntry("LiO2", -15.0, entry_id="e4").as_dict()]
    with open(os.path.join(d, "entries.json"), "w") as fh:
        json.dump(ents, fh)
    with open(os.path.join(d, "hyp.txt"), "w") as fh:
        fh.write("e4\n")
    p = subprocess.run(CLI + ["convert-entries", "--in", os.path.join(d, "entries.json"), "--out",
                              os.path.join(d, "conv.csv"), "--hypothetical-ids", os.path.join(d, "hyp.txt")],
                       capture_output=True, text=True)
    rows = read_csv(os.path.join(d, "conv.csv")) if p.returncode == 0 else []
    check("convert-entries", p.returncode == 0 and len(rows) == 4 and sum(r["is_hypothetical"] == "1" for r in rows) == 1,
          p.stderr[-200:] if p.returncode else f"{len(rows)} rows")
    snap = read_snapshot(os.path.join(d, "snapshots", "srcb.csv"))
    one = compute_dhd(snap, lambda r: r["is_hypothetical"], workers=1)
    two = compute_dhd(snap, lambda r: r["is_hypothetical"], workers=2)
    same = all((math.isnan(one[k][0]) and math.isnan(two[k][0])) or abs(one[k][0] - two[k][0]) < 1e-12 for k in one)
    check("parallel hull equals serial hull", same and len(one) == len(two), f"{len(one)} entries")
    cfgp = os.path.join(d, "config.json")
    cfg = json.load(open(cfgp))
    raw = [ComputedEntry(r["formula"], r_e * 1, entry_id=r["entry_id"]).as_dict()
           for r, r_e in ((x, float(x["energy_per_atom"])) for x in read_csv(os.path.join(d, "snapshots", "srca_r2.csv")))]
    with open(os.path.join(d, "raw_r2.json"), "w") as fh:
        json.dump(raw, fh)
    cfg["relabeling"] = [{"name": "smoke", "release_r": "srca_r1", "release_r1": "srca_r2", "raw_entries_r1": "raw_r2.json",
                          "scheme_r": "MaterialsProjectCompatibility", "scheme_r1": "MaterialsProject2020Compatibility",
                          "corrected_species": ["O", "F"]}]
    cfg["cutoff_date"] = "2026-12-31"; cfg["allow_future_cutoff_for_testing"] = True
    json.dump(cfg, open(cfgp, "w"))
    p = run(d, "relabel-pop", ok=False)
    rows = read_csv(os.path.join(d, "work", "cohorts", "relabel__smoke.csv")) if p.returncode == 0 else []
    check("relabel-pop runs (entries without calculation parameters are unclassifiable)",
          p.returncode == 0 and rows and all(r["group"] in {"M", "flip_unclassifiable"} for r in rows),
          (p.stderr[-300:] if p.returncode else f"{len(rows)} rows"))


def oqmd_converter_test(d):
    """convert-oqmd on a miniature export in the exact OQMD column layout, then hull and x0."""
    H = ["entry_id", "composition_id", "path", "label", "duplicate_of_id", "prototype_id", "ntypes", "natoms",
         "entry_delta_e", "entry_stability", "fe_id", "fit_id", "fe_delta_e", "fe_stability"]
    P = "/home/oqmd/libraries/"
    rows = [
        # elements (formation energy near 0 against the fitted reference)
        [1, "Li1", P + "prototypes/elements/A2/Li", "A2", "NULL", "A2", 1, 2, "NULL", "NULL", 11, "standard", 0.004, 0.004],
        [2, "O1", P + "prototypes/elements/A1/O", "A1", "NULL", "A1", 1, 2, "NULL", "NULL", 12, "standard", 0.0, 0.0],
        [3, "Fe1", P + "prototypes/elements/A2/Fe", "A2", "NULL", "A2", 1, 2, "NULL", "NULL", 13, "standard", 0.0, 0.0],
        # ICSD-origin Li2O
        [4, "Li2 O1", P + "icsd/101", "icsd-101", 4, "C1", 2, 3, "NULL", "NULL", 14, "standard", -2.0, -0.05],
        # hypothetical decoration at the ICSD composition (source-known)
        [5, "Li2 O1", P + "prototypes/binaries/B1/Li_O", "NULL", "NULL", "B1", 2, 3, "NULL", "NULL", 15, "standard", -1.9, 0.10],
        # hypothetical, stable, new (half-heusler path, NULL prototype as in v1.0), two fe rows, latest wins
        [6, "Fe1 Li1 O1", P + "prototypes/ternaries/half_heusler/Fe/Li_O", "NULL", "NULL", "NULL", 3, 3, "NULL", "NULL", 16, "standard", -1.0, 0.02],
        [6, "Fe1 Li1 O1", P + "prototypes/ternaries/half_heusler/Fe/Li_O", "NULL", "NULL", "NULL", 3, 3, "NULL", "NULL", 26, "standard", -1.6, -0.03],
        # duplicate of entry 6 (dropped)
        [7, "Fe1 Li1 O1", P + "prototypes/ternaries/heusler/Fe/Li_O", "NULL", 6, "NULL", 3, 4, "NULL", "NULL", 17, "standard", -1.6, -0.03],
        # hypothetical binary near the hull, displayed label unknown
        [8, "Fe1 O1", P + "prototypes/binaries/B1/Fe_O", "NULL", "NULL", "B1", 2, 2, "NULL", "NULL", 18, "standard", -1.3, "NULL"],
        # no formation energy (dropped), test path (kept in hull, not hypothetical)
        [9, "Fe2 O3", P + "prototypes/binaries/D5_1/Fe_O", "NULL", "NULL", "D5_1", 2, 5, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL"],
        [10, "Fe1 Li2", "/home/oqmd/test/tmp/x", "NULL", "NULL", "NULL", 2, 3, "NULL", "NULL", 20, "standard", 0.2, 0.2],
        # ICSD entry with no formation energy: its composition is still known to the source
        [11, "Fe3 O4", P + "icsd/202", "icsd-202", 11, "H1_1", 2, 14, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL"],
        [12, "Fe3 O4", P + "prototypes/binaries/H1_1/Fe_O", "NULL", "NULL", "H1_1", 2, 7, "NULL", "NULL", 22, "standard", -1.0, 0.5],
        # hypothetical with two rows where the LATEST row is the higher one: the lowest row must be kept
        [13, "Li1 O2", P + "prototypes/binaries/C2/Li_O", "NULL", "NULL", "C2", 2, 3, "NULL", "NULL", 23, "standard", -1.5, -0.01],
        [13, "Li1 O2", P + "prototypes/binaries/C2/Li_O", "NULL", "NULL", "C2", 2, 3, "NULL", "NULL", 33, "standard", -1.2, 0.29],
    ]
    os.makedirs(os.path.join(d, "oq"), exist_ok=True)
    with open(os.path.join(d, "oq", "e.tsv"), "w") as fh:
        fh.write("\t".join(H) + "\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    with open(os.path.join(d, "oq", "s.tsv"), "w") as fh:
        fh.write("entry_id\tstructure_id\tlabel\tspacegroup_id\n6\t1\tinput\t216\n6\t2\tstandard\t225\n")
    p = subprocess.run(CLI + ["convert-oqmd", "--entries", os.path.join(d, "oq", "e.tsv"), "--structures",
                              os.path.join(d, "oq", "s.tsv"), "--out", os.path.join(d, "oq", "snap.csv")],
                       capture_output=True, text=True)
    rep = read_json(os.path.join(d, "oq", "snap_conversion_report.json")) if p.returncode == 0 else {}
    snap = {r["entry_id"]: r for r in read_csv(os.path.join(d, "oq", "snap.csv"))} if p.returncode == 0 else {}
    ok = (p.returncode == 0 and rep["dropped_duplicates"] == 1 and rep["dropped_no_formation_energy"] == 2
          and rep["entries_with_conflicting_delta_e_gt_1meV"] == 2 and snap["oqmd-6"]["energy_per_atom"] == "-1.6"
          and snap["oqmd-13"]["energy_per_atom"] == "-1.5" and snap["oqmd-13"]["displayed_stable"] == "1"
          and rep["entries_where_lowest_row_is_not_latest_row"] == 1
          and snap["oqmd-12"]["source_known_composition"] == "1" and rep["icsd_compositions_known_to_source"] == 2
          and snap["oqmd-6"]["spacegroup_number"] == "225" and "heusler" in snap["oqmd-6"]["prototype"]
          and snap["oqmd-5"]["source_known_composition"] == "1" and snap["oqmd-4"]["is_hypothetical"] == "0"
          and snap["oqmd-10"]["is_hypothetical"] == "0" and snap["oqmd-8"]["displayed_stable"] == ""
          and rep["elemental_references_added"] == 3 and rep["displayed_stable_new_heusler_type_compositions"] == 1
          and rep["displayed_stable_new_compositions"] == 2)
    check("convert-oqmd applies the diagnostics-based rules", ok, p.stderr[-300:] if p.returncode else f"{rep}")
    cfg = {"workdir": "work_oq", "cutoff_date": "2026-12-11", "allow_future_cutoff_for_testing": True,
           "releases": [{"source": "OQMD", "release_id": "v1_0", "release_date": "2014-04", "path": "oq/snap.csv",
                         "granularity": "composition", "tau_ev_per_atom": 0.0,
                         "expected_counts": {"stable_new_compositions": 2, "stable_new_heusler_type_compositions": 1}}],
           "threshold": {"sources": ["OQMD"], "windows_mev": [10, 25, 50], "horizon_years": 5,
                         "assumed_outcome_rate": 0.05, "effect_of_interest": 0.03},
           "seeds": {"expert_sample": 1, "two_phase": 2, "residual_audit": 3, "rerun": 4}}
    with open(os.path.join(d, "oq_config.json"), "w") as fh:
        json.dump(cfg, fh)
    c2 = CLI + ["--config", os.path.join(d, "oq_config.json")]
    outs = [subprocess.run(c2 + [step], capture_output=True, text=True) for step in ("check", "hull", "x0", "cohorts")]
    x0 = read_json(os.path.join(d, "work_oq", "checks", "x0.json")) if all(o.returncode == 0 for o in outs) else {}
    rec = x0.get("v1_0", {})
    coh = read_csv(os.path.join(d, "work_oq", "cohorts", "oqmd.csv")) if rec else []
    thr = {r["reduced_formula"]: r for r in read_csv(os.path.join(d, "work_oq", "cohorts", "threshold__oqmd.csv"))} if rec else {}
    check("threshold running variable is the displayed decomposition energy",
          "LiO2" in thr and abs(float(thr["LiO2"]["x_mev"]) + 10.0) < 1e-6 and thr["LiO2"]["x_recomputed_mev"] != ""
          and rec.get("new_composition_overlap", {}).get("jaccard") is not None,
          f"{ {k: (v['x_mev'], v['x_recomputed_mev'], v['side']) for k, v in thr.items()} }; overlap {rec.get('new_composition_overlap')}")
    check("OQMD snapshot: hull, x0 and cohort (source-known compositions excluded)",
          rec.get("within_tolerance") is True and sorted(c["reduced_formula"] for c in coh) == ["LiFeO", "LiO2"]
          and rec["displayed"]["stable_new_compositions"] == 2
          and all(c["entry_date"] == "2014-04-15" for c in coh),
          f"x0 {rec.get('displayed')} / {rec.get('recomputed')}; cohort {[c['reduced_formula'] for c in coh]}"
          + "".join(o.stderr[-200:] for o in outs if o.returncode))


def mp_converter_test(d):
    """convert-mp on a miniature build in the 2022-10-28 layout, then the premise check
    against the miniature OQMD snapshot from oqmd_converter_test."""
    import gzip
    B = os.path.join(d, "mpb")
    def rec(mid, comp, fe, stable, dhd, tt="GGA_GGA+U"):
        return {"material_id": mid, "thermo_type": tt, "composition_reduced": comp, "formation_energy_per_atom": fe,
                "is_stable": stable, "decomposition_enthalpy": dhd, "energy_above_hull": max(dhd, 0.0),
                "deprecated": False, "symmetry": {"number": 225}, "entries": {"GGA": {"big": "x" * 50}}}
    thermo = [rec("mp-1", {"Li": 2.0, "O": 1.0}, -2.1, True, -0.30),
              rec("mp-2", {"Li": 1.0, "Fe": 1.0, "O": 1.0}, -1.55, False, 0.02),
              rec("mp-3", {"Li": 1.0, "O": 2.0}, -1.4, False, 0.05),
              rec("mp-1", {"Li": 2.0, "O": 1.0}, -2.0, True, -0.2, tt="GGA_GGA+U_R2SCAN"),
              rec("mp-4", {"Fe": 1.0}, 0.0, True, 0.0)]
    prov = [{"material_id": "mp-1", "theoretical": False, "created_at": {"$date": "2022-10-20T00:00:00Z"},
             "database_IDs": {"icsd": ["icsd-1"]}},
            {"material_id": "mp-2", "theoretical": True, "created_at": {"$date": "2022-10-20T00:00:00Z"}, "database_IDs": {}},
            {"material_id": "mp-3", "theoretical": True, "created_at": {"$date": "2019-01-01T00:00:00Z"}, "database_IDs": {}}]
    os.makedirs(os.path.join(B, "thermo", "thermo_type=GGA_GGA+U", "nelements=2"), exist_ok=True)
    os.makedirs(os.path.join(B, "provenance", "nelements=2"), exist_ok=True)
    with gzip.open(os.path.join(B, "thermo", "thermo_type=GGA_GGA+U", "nelements=2", "a.jsonl.gz"), "wt") as fh:
        fh.write("\n".join(json.dumps(r) for r in thermo) + "\n")
    with gzip.open(os.path.join(B, "provenance", "nelements=2", "p.jsonl.gz"), "wt") as fh:
        fh.write("\n".join(json.dumps(r) for r in prov) + "\n")
    out = os.path.join(d, "oq", "mp.csv")
    p = subprocess.run(CLI + ["convert-mp", "--build-dir", B, "--out", out], capture_output=True, text=True)
    rep = read_json(out.replace(".csv", "_conversion_report.json")) if p.returncode == 0 else {}
    snap = {r["entry_id"]: r for r in read_csv(out)} if p.returncode == 0 else {}
    ok = (p.returncode == 0 and rep["chosen_materials"] == 4 and snap["mp-1"]["energy_per_atom"] == "-2.1"
          and snap["mp-1"]["origin_icsd_id"] == "icsd" and snap["mp-2"]["is_hypothetical"] == "1"
          and snap["mp-2"]["displayed_dhd"] == "0.02" and rep["provenance_created_at_years_theoretical"] == {"2019": 1, "2022": 1}
          and snap["mp-4"]["is_hypothetical"] == "0")
    check("convert-mp chooses records by their thermo_type field and reads provenance", ok,
          p.stderr[-300:] if p.returncode else f"{ {k: rep[k] for k in ('chosen_materials', 'hypothetical', 'thermo_records_by_folder_and_field')} }")
    cfg = json.load(open(os.path.join(d, "oq_config.json")))
    cfg["releases"].append({"source": "Materials Project", "release_id": "mp_2022", "release_date": "2022-10-28",
                            "path": "oq/mp.csv", "granularity": "structure", "tau_ev_per_atom": 0.0,
                            "use_for_cohorts": False})
    cfg["premise"] = {"pairs": [{"a": "v1_0", "b": "mp_2022", "match": "composition"}]}
    json.dump(cfg, open(os.path.join(d, "oq_config.json"), "w"))
    c2 = CLI + ["--config", os.path.join(d, "oq_config.json")]
    outs = [subprocess.run(c2 + [step], capture_output=True, text=True) for step in ("cohorts", "premise")]
    prem = read_json(os.path.join(d, "work_oq", "checks", "premise.json")) if all(o.returncode == 0 for o in outs) else {}
    r = prem.get("v1_0 vs mp_2022", {})
    check("premise check from displayed values; premise-only release kept out of cohorts",
          r.get("matched") == 3 and not os.path.exists(os.path.join(d, "work_oq", "cohorts", "materials_project.csv")),
          f"{r}" + "".join(o.stderr[-300:] for o in outs if o.returncode))


def cod_converter_test(d):
    """convert-cod on rows in the COD `data` layout; then dbmatch reads its output."""
    from ncr.cod import FIELDS
    rows = [
        {"file": "1000001", "formula": "- Fe3 O4 -", "year": "1986", "method": "", "status": "", "flags": "has coordinates",
         "sgNumber": "227", "cellpressure": "", "diffrpressure": ""},
        {"file": "1000002", "formula": "- Li Fe O -", "year": "2019", "method": "single crystal", "status": "",
         "flags": "has coordinates,has disorder", "sgNumber": "225", "cellpressure": "2000000", "diffrpressure": ""},
        {"file": "1000003", "formula": "- Li O2 -", "year": "2020", "method": "theoretical", "status": ""},
        {"file": "1000004", "formula": "- Fe O -", "year": "2001", "method": "", "status": "retracted"},
        {"file": "1000005", "formula": "", "calcformula": "- Li2 O -", "year": "", "method": "", "status": ""},
        {"file": "1000006", "formula": "- C6 H6 -", "year": "2005", "method": "", "status": ""},
    ]
    p = os.path.join(d, "cod.tsv")
    with open(p, "w") as fh:
        fh.write("\t".join(FIELDS) + "\n")
        for r in rows:
            fh.write("\t".join(r.get(k, "") for k in FIELDS) + "\n")
    out = os.path.join(d, "cod_exp.csv")
    q = subprocess.run(CLI + ["convert-cod", "--in", p, "--out", out], capture_output=True, text=True)
    rep = read_json(out.replace(".csv", "_conversion_report.json")) if q.returncode == 0 else {}
    ex = {r["exp_id"]: r for r in read_csv(out)} if q.returncode == 0 else {}
    ok = (q.returncode == 0 and set(ex) == {"cod-1000001", "cod-1000002", "cod-1000006"}
          and rep.get("excluded_theoretical") == 1 and rep.get("excluded_retracted") == 1
          and rep.get("excluded_no_year") == 1 and ex["cod-1000002"]["high_pressure"] == "1"
          and ex["cod-1000002"]["disordered"] == "1" and ex["cod-1000001"]["first_report_date"] == "1986"
          and rep.get("contains_C_and_H") == 1)
    check("convert-cod applies the COD rules", ok, q.stderr[-300:] if q.returncode else f"{rep}")
    from ncr.cohorts import read_experiments
    exps = read_experiments(out)
    check("COD table is readable by dbmatch", len(exps) == 3 and str(exps[0]["date"]) == "1986-07-01",
          f"{[(e['exp_id'], str(e['date'])) for e in exps]}")
    cfg = json.load(open(os.path.join(d, "oq_config.json")))
    cfg["experiments"] = {"path": "cod_exp.csv"}
    json.dump(cfg, open(os.path.join(d, "oq_config.json"), "w"))
    q = subprocess.run(CLI + ["--config", os.path.join(d, "oq_config.json"), "coverage-known"], capture_output=True, text=True)
    cov = read_json(os.path.join(d, "work_oq", "checks", "coverage_known.json")) if q.returncode == 0 else {}
    check("coverage-known reads only source-known compounds", cov.get("v1_0", {}).get("icsd_origin_compositions") == 1
          and cov.get("mp_2022", {}).get("icsd_origin_compositions") == 1,
          q.stderr[-300:] if q.returncode else f"{ {k: v for k, v in cov.items() if k != 'experimental_reports_compositions'} }")


def census_test(d, truth):
    """Annotation-free design: two dated sources, dbmatch with per-source dates, census
    estimates, capture-recapture; plus the two ICSD-derived converters."""
    import gzip
    exps = read_csv(os.path.join(d, "experiments.csv"))
    rng = np.random.default_rng(11)
    a_rows, b_rows = [], []
    for r in exps:
        u = rng.random()
        if u < 0.7:
            a_rows.append(dict(r, database="COD", exp_id="c" + r["exp_id"]))
        if u > 0.4:
            b_rows.append(dict(r, database="MP-ICSD", exp_id="m" + r["exp_id"]))
    write_csv(os.path.join(d, "srcA.csv"), a_rows); write_csv(os.path.join(d, "srcB.csv"), b_rows)
    q = subprocess.run(CLI + ["merge-experiments", "--in", os.path.join(d, "srcA.csv"), os.path.join(d, "srcB.csv"),
                              "--out", os.path.join(d, "merged.csv")], capture_output=True, text=True)
    cfg = json.load(open(os.path.join(d, "config.json")))
    cfg["workdir"] = "work_census"; cfg["experiments"] = {"path": "merged.csv"}
    json.dump(cfg, open(os.path.join(d, "census_config.json"), "w"))
    c2 = CLI + ["--config", os.path.join(d, "census_config.json")]
    outs = [subprocess.run(c2 + st, capture_output=True, text=True)
            for st in (["check"], ["hull"], ["cohorts"], ["dbmatch", "--cohort", "all"], ["estimate-census"])]
    ok_run = q.returncode == 0 and all(o.returncode == 0 for o in outs)
    res = read_json(os.path.join(d, "work_census", "results", "census_results.json")) if ok_run else {}
    db = read_csv(os.path.join(d, "work_census", "dbmatch", "srcb.csv")) if ok_run else []
    rows = read_csv(os.path.join(d, "work_census", "cohorts", "srcb.csv")) if ok_run else []
    e = {r["pred_id"]: decimal_year(parse_date(r["entry_date"])[0]) for r in rows}
    direct = sum(1 for x in db if x["db_first_T3"] and decimal_year(parse_date(x["db_first_T3"])[0]) <= e[x["pred_id"]])
    pe = res.get("base", {}).get("srcb", {}).get("pre_existence_all", {})
    check("census pre-existence equals the direct count from matched records",
          ok_run and pe.get("k") == direct and pe.get("n") == len(rows),
          "".join(o.stderr[-200:] for o in outs if o.returncode) or f"{pe.get('k')} vs {direct} of {len(rows)}")
    srcs_ok = ok_run and "db_first_T3__COD" in db[0] and "db_first_T3__MP-ICSD" in db[0]
    cr = res.get("capture_recapture", {}).get("srcb", {})
    check("per-source dates and capture-recapture", srcs_ok and cr.get("chapman_total", 0) >= cr.get("union", 1),
          f"{ {k: (round(v, 3) if isinstance(v, float) else v) for k, v in cr.items() if k != 'assumption'} }")
    q3 = subprocess.run(c2 + ["robustness"], capture_output=True, text=True)
    rob = read_json(os.path.join(d, "work_census", "robustness_for_claude.json")) if q3.returncode == 0 else {}
    need = {"composition_tolerance_halved", "composition_tolerance_doubled", "without_high_pressure_reports",
            "year_only_dates_1_january", "year_only_dates_31_december"}
    rr = rob.get("reanalysis", {}).get("threshold__srcb", {})
    check("robustness command runs every specified variant",
          q3.returncode == 0 and need <= set(rob.get("rematch", {})) and "local_linear_25meV" in rr
          and "one_entry_per_chemical_system" in rr and "blocks_without_number_of_elements" in rr,
          q3.stderr[-400:] if q3.returncode else f"rematch {sorted(rob['rematch'])}; reanalysis {sorted(rr)}")
    if rob:
        p_half = rob["rematch"]["composition_tolerance_halved"]["srcb"]["srcb::srcb"]["pre_existence"]["k"]
        p_dbl = rob["rematch"]["composition_tolerance_doubled"]["srcb"]["srcb::srcb"]["pre_existence"]["k"]
        p0 = rob["primary"]["srcb"]["srcb::srcb"]["pre_existence"]["k"]
        check("tolerance variants are nested (halved <= primary <= doubled)", p_half <= p0 <= p_dbl, f"{p_half} <= {p0} <= {p_dbl}")
    check("census rates per entry release", any(k.startswith("srca::") for k in res.get("base", {})),
          f"{[k for k in res.get('base', {})]}")
    th = res.get("threshold", {}).get("threshold__srcb", {}).get("windows", {}).get("25", {})
    check("threshold effect from a complete census keeps the assignment variance",
          "effect" in th and th["effect"]["se"] > 0 and abs(th["effect"]["se_certification_only"]) < 1e-12
          and "attention_check" not in th, f"{th.get('effect')}")
    # converters for the ICSD-derived sources
    with open(os.path.join(d, "oqi.tsv"), "w") as fh:
        fh.write("entry_id\tcomposition_id\tyear\n1\tFe3 O4\t2010\n2\tLi2 O1\tNULL\n3\tNa1 Cl1\t1995\n")
    q1 = subprocess.run(CLI + ["convert-oqmd-icsd", "--in", os.path.join(d, "oqi.tsv"), "--out", os.path.join(d, "oqi.csv")],
                        capture_output=True, text=True)
    r1 = read_json(os.path.join(d, "oqi_conversion_report.json")) if q1.returncode == 0 else {}
    B = os.path.join(d, "mpb2")
    os.makedirs(os.path.join(B, "provenance", "x"), exist_ok=True); os.makedirs(os.path.join(B, "thermo", "x"), exist_ok=True)
    prov = [{"material_id": "mp-1", "theoretical": False, "database_IDs": {"icsd": ["icsd-1"]},
             "references": ["@article{a, year = {2005}}", "@article{b, year = {1999}}"]},
            {"material_id": "mp-2", "theoretical": False, "database_IDs": {}, "references": []},
            {"material_id": "mp-3", "theoretical": True, "database_IDs": {}, "references": []}]
    therm = [{"material_id": "mp-1", "composition_reduced": {"Fe": 3.0, "O": 4.0}},
             {"material_id": "mp-2", "composition_reduced": {"Li": 2.0, "O": 1.0}}]
    with gzip.open(os.path.join(B, "provenance", "x", "p.jsonl.gz"), "wt") as fh:
        fh.write("\n".join(json.dumps(x) for x in prov) + "\n")
    with gzip.open(os.path.join(B, "thermo", "x", "t.jsonl.gz"), "wt") as fh:
        fh.write("\n".join(json.dumps(x) for x in therm) + "\n")
    q2 = subprocess.run(CLI + ["convert-mp-icsd", "--build-dir", B, "--out", os.path.join(d, "mpi.csv")],
                        capture_output=True, text=True)
    r2 = read_json(os.path.join(d, "mpi_conversion_report.json")) if q2.returncode == 0 else {}
    mpi = read_csv(os.path.join(d, "mpi.csv")) if q2.returncode == 0 else []
    check("ICSD-derived converters (OQMD-ICSD, MP-ICSD) apply their rules",
          r1.get("rows_written") == 2 and r1.get("year_max") == 2010 and r2.get("rows_written") == 1
          and r2.get("excluded_no_icsd_identifier") == 1 and mpi and mpi[0]["first_report_date"] == "1999",
          (q1.stderr + q2.stderr)[-300:] or f"{r1} | {r2}")


def density_test_large_n():
    """The density (binomial) test must work for thousands of entries and match small exact cases."""
    from ncr.estimate import _binom_two_sided
    p_mid, p_far, p_small = _binom_two_sided(2000, 4000), _binom_two_sided(1800, 4000), _binom_two_sided(0, 10)
    check("density test is exact for small n and stable for large n",
          abs(p_mid - 1.0) < 1e-9 and 0 < p_far < 1e-8 and abs(p_small - 2 / 1024) < 1e-12,
          f"p(2000/4000)={p_mid:.6f}, p(1800/4000)={p_far:.3e}, p(0/10)={p_small:.6f}")


def polymorph_test():
    """A higher-energy polymorph of a stable composition must get its energy above the hull."""
    from ncr.snapshots import compute_dhd
    rows = [{"entry_id": "li", "formula": "Li", "energy_per_atom": 0.0, "is_hypothetical": False},
            {"entry_id": "o", "formula": "O", "energy_per_atom": 0.0, "is_hypothetical": False},
            {"entry_id": "gs", "formula": "Li2O", "energy_per_atom": -2.0, "is_hypothetical": True},
            {"entry_id": "poly", "formula": "Li2O", "energy_per_atom": -1.9, "is_hypothetical": True}]
    r = compute_dhd(rows, lambda x: x["is_hypothetical"])
    check("higher-energy polymorph scored above the hull, ground state below",
          abs(r["poly"][0] - 0.1) < 1e-9 and r["gs"][0] < 0, f"gs {r['gs'][0]:.3f}, poly {r['poly'][0]:.3f}")


def two_phase_ratio_coverage(R=400):
    """Coverage of the two-phase ratio interval (delta method) by simulation."""
    import nc_realization_metrics as NCM
    from ncr.estimate import two_phase_ratio_adjust
    rng = np.random.default_rng(3)
    Npop, n1 = 20000, 2000
    P = rng.random(Npop) < 0.15
    Y = (~P) & (rng.random(Npop) < 0.12)
    trueR = Y.sum() / (~P).sum()
    cover = 0
    for _ in range(R):
        s1 = rng.choice(Npop, n1, replace=False)
        y, dd = Y[s1].astype(float), (~P[s1]).astype(float)
        Ly = np.where(rng.random(n1) < 0.9, y, 1 - y)
        Ld = np.where(rng.random(n1) < 0.95, dd, 1 - dd)
        Ly = Ly * Ld
        strata = np.where(Ly == 1, "pos", "neg").astype(object)
        idx = np.concatenate([rng.choice(np.flatnonzero(strata == k), 150, replace=False) for k in ("pos", "neg")])

        class B:
            N = n1

            def cd(self, c, L, H):
                return NCM.census_difference(c, L, strata, idx, H)
        ones = np.ones(n1)
        rr = NCM.census_ratio(ones, ones, Ly, Ld, strata, idx, y[idx], dd[idx])
        adj = two_phase_ratio_adjust(rr, B(), Ly, y[idx], Ld, dd[idx], Npop, 0.05)
        cover += adj["ci_low"] <= trueR <= adj["ci_high"]
    cov = cover / R
    check("two-phase ratio interval coverage", 0.92 <= cov <= 0.985, f"{cov:.3f} over {R} simulations (target 0.95)")


if __name__ == "__main__":
    main()
