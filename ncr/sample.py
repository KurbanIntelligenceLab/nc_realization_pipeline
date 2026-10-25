"""Expert sampling (after freezing), blinded worksheets, adjudication and the
residual-miss audit sample."""
from __future__ import annotations

import datetime as dt
import os
import sys

import numpy as np

from .util import date_str, decimal_year, log, parse_date, read_csv, sha256_file, write_csv, write_json

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)
import nc_realization_metrics as NCM                 # noqa: E402

EXPERT_FIELDS = ["first_report_date", "first_report_received_date", "first_report_reference", "tier",
                 "high_pressure", "thin_film_only", "structure_decision_without_refined_structure",
                 "first_mention_date", "minutes_spent", "sources_searched", "notes"]
BLIND_FIELDS = ["sample_id", "formula", "elements", "granularity", "structure_file",
                "predicted_spacegroup", "predicted_prototype", "entry_date"]


def strata_for(labs: list[dict]) -> np.ndarray:
    """Strata fixed before the draw: cohort x model label, small strata merged (Table S1).
    The estimation step recomputes them with this same function."""
    strata = np.array([f"{x['cohort']}|{x['model_label']}" for x in labs], dtype=object)
    return NCM.merge_small_strata(strata, 2)


def _dirs(cfg, labels_file):
    name = os.path.splitext(os.path.basename(labels_file))[0]
    d = os.path.join(cfg["workdir"], "experts", name)
    os.makedirs(d, exist_ok=True)
    return name, d


def step_sample(cfg: dict, cohort_file: str, labels_file: str, manifest, n_total: int | None = None) -> str:
    frozen = labels_file + ".sha256"
    if not os.path.exists(frozen):
        raise RuntimeError("labels are not frozen; run 'freeze' before drawing the expert sample")
    h = sha256_file(labels_file)
    if open(frozen, encoding="utf-8").read().strip() != h:
        raise RuntimeError("labels changed after freezing; the sample cannot be drawn")
    name, d = _dirs(cfg, labels_file)
    if os.path.exists(os.path.join(d, "sample_key.csv")):
        raise RuntimeError(f"a sample already exists for {name}; delete {d} only if no expert has started")
    sc = cfg.get("sampling", {})
    n_total = int(n_total or sc.get("n_per_cohort", 400))
    min_per = int(sc.get("min_per_stratum", 30))
    rates = sc.get("assumed_error_rates", {"report": 0.2, "none": 0.05, "abstain": 0.3})
    seed = int(cfg["seeds"]["expert_sample"])
    rows = {r["pred_id"]: r for r in read_csv(cohort_file)}
    labs = read_csv(labels_file)
    strata = strata_for(labs)
    names, counts = np.unique(strata, return_counts=True)
    S = np.array([np.sqrt(rates.get(n.split("|")[-1], 0.2) * (1 - rates.get(n.split("|")[-1], 0.2))) for n in names])
    n_h = NCM.neyman_allocation(counts, S, n_total, min_per)
    n_h = np.maximum(np.minimum(n_h, counts), np.minimum(2, counts))
    rng = np.random.default_rng(seed)
    idx = NCM.stratified_srswor(strata, {k: int(v) for k, v in zip(names, n_h)}, rng)
    ids, seen = [], set()
    while len(ids) < idx.size:
        v = int(rng.integers(0, 16 ** 8))
        if v not in seen:
            seen.add(v)
            ids.append(v)
    key, sheet = [], []
    for j, i in enumerate(idx):
        lab = labs[i]
        r = rows[lab["pred_id"]]
        sid = f"S{ids[j]:08x}"
        key.append({"sample_id": sid, "pred_id": lab["pred_id"], "cohort": lab["cohort"], "stratum": strata[i],
                    "row_index": int(i)})
        sheet.append({"sample_id": sid, "formula": r["formula"], "elements": r["elements"],
                      "granularity": r.get("granularity", ""), "structure_file": r.get("structure_file", ""),
                      "predicted_spacegroup": r.get("spacegroup_number", ""),
                      "predicted_prototype": r.get("prototype", ""), "entry_date": r["entry_date"],
                      **{f: "" for f in EXPERT_FIELDS}})
    write_csv(os.path.join(d, "sample_key.csv"), key)
    for who in ("A", "B"):
        order = rng.permutation(len(sheet))
        write_csv(os.path.join(d, f"expert_{who}.csv"), [sheet[k] for k in order], BLIND_FIELDS + EXPERT_FIELDS)
    alloc = {str(k): {"N": int(c), "n": int(n)} for k, c, n in zip(names, counts, n_h)}
    write_json(os.path.join(d, "allocation.json"), {"labels_sha256": h, "seed": seed, "allocation": alloc})
    manifest.add("expert_sample", labels_file=os.path.relpath(labels_file, cfg["workdir"]), labels_sha256=h,
                 seed=seed, n=int(idx.size), allocation=alloc)
    log(f"sample: {name}: {idx.size} entries; give expert_A.csv and expert_B.csv to the two experts "
        f"(never sample_key.csv)")
    return d


# ----------------------------------------------------------------------------
def _cutpoints(entry: dt.date, horizons, extra=()):
    e = decimal_year(entry)
    return sorted({e, *[e + t for t in horizons], *extra})


def _bucket(date_s: str, cuts) -> int:
    d, _ = parse_date(date_s)
    if d is None:
        return -1
    y = decimal_year(d)
    return int(sum(1 for c in cuts if y > c))


def _qualifies(tier: str, gran: str) -> bool:
    t = (tier or "").strip().upper()
    return t in ({"T1", "T2"} if gran == "structure" else {"T1", "T2", "T3"})


def _expert_view(row: dict, gran: str) -> tuple[str, str]:
    """(first qualifying report date or '', first mention date or '')."""
    rep = row["first_report_date"] if _qualifies(row.get("tier", ""), gran) else ""
    return rep, row.get("first_mention_date", "")


def step_adjudicate(cfg: dict, cohort_file: str, labels_file: str, finalize: bool = False) -> str:
    name, d = _dirs(cfg, labels_file)
    key = {k["sample_id"]: k for k in read_csv(os.path.join(d, "sample_key.csv"))}
    rows = {r["pred_id"]: r for r in read_csv(cohort_file)}
    A = {x["sample_id"]: x for x in read_csv(os.path.join(d, "expert_A.csv"))}
    B = {x["sample_id"]: x for x in read_csv(os.path.join(d, "expert_B.csv"))}
    if set(A) != set(key) or set(B) != set(key):
        raise ValueError("expert sheets do not contain exactly the sampled entries")
    horizons = sorted({*cfg.get("horizons_years", [3, 5, 10]), *cfg.get("ml_horizons_years", [2]),
                       float(cfg.get("threshold", {}).get("horizon_years", 5))})
    adj_path = os.path.join(d, "adjudication.csv")
    if finalize:
        adj = {x["sample_id"]: x for x in read_csv(adj_path)}
    out, kap_a, kap_b, kap3_a, kap3_b, todo = [], [], [], [], [], []
    for sid, k in key.items():
        r = rows[k["pred_id"]]
        gran = r.get("granularity", "structure")
        entry, _ = parse_date(r["entry_date"])
        extra = []
        if r.get("release_r_date"):
            extra.append(decimal_year(parse_date(r["release_r_date"])[0]))
        cuts = _cutpoints(entry, horizons, extra)
        ra, ma = _expert_view(A[sid], gran)
        rb, mb = _expert_view(B[sid], gran)
        for v, who in ((A[sid], "A"), (B[sid], "B")):
            if v.get("tier", "").strip().upper() not in {"T1", "T2", "T3", "NONE"}:
                raise ValueError(f"expert {who}, {sid}: tier must be T1, T2, T3 or none")
            if v.get("tier", "").strip().upper() != "NONE" and not v.get("first_report_date"):
                raise ValueError(f"expert {who}, {sid}: a tier without a first_report_date")
        kap_a.append(int(bool(ra)))
        kap_b.append(int(bool(rb)))
        cat = lambda s: "none" if not s else ("pre" if decimal_year(parse_date(s)[0]) <= decimal_year(entry) else "later")
        kap3_a.append(cat(ra))
        kap3_b.append(cat(rb))
        agree = (_bucket(ra, cuts) == _bucket(rb, cuts)) and (_bucket(ma, cuts) == _bucket(mb, cuts))
        base = {"sample_id": sid, "pred_id": k["pred_id"], "cohort": k["cohort"]}
        if agree:
            pick = lambda x, y: min([s for s in (x, y) if s], default="")
            out.append({**base, "final_first_report_date": pick(ra, rb),
                        "final_first_received_date": pick(A[sid].get("first_report_received_date", "") if ra else "",
                                                          B[sid].get("first_report_received_date", "") if rb else ""),
                        "final_first_mention_date": pick(ma, mb), "adjudicated_by": "agreement"})
        elif finalize:
            x = adj[sid]
            if x.get("final_first_report_date", "") == "" and x.get("final_has_report", "").lower() not in {"no", "n", "0"}:
                raise ValueError(f"{sid}: adjudicator must give final_first_report_date or final_has_report=no")
            out.append({**base, "final_first_report_date": x.get("final_first_report_date", ""),
                        "final_first_received_date": x.get("final_first_received_date", ""),
                        "final_first_mention_date": x.get("final_first_mention_date", ""),
                        "adjudicated_by": "third_expert"})
        else:
            todo.append({**base, "formula": r["formula"], "entry_date": r["entry_date"],
                         "A_report_date": ra, "A_reference": A[sid].get("first_report_reference", ""),
                         "A_tier": A[sid].get("tier", ""), "A_mention": ma,
                         "B_report_date": rb, "B_reference": B[sid].get("first_report_reference", ""),
                         "B_tier": B[sid].get("tier", ""), "B_mention": mb,
                         "final_has_report": "", "final_first_report_date": "", "final_first_received_date": "",
                         "final_first_mention_date": "", "adjudicator_notes": ""})
    kappa = {"has_report": NCM.cohen_kappa(np.array(kap_a), np.array(kap_b)),
             "pre_later_none": NCM.cohen_kappa(np.array(kap3_a), np.array(kap3_b)), "n": len(key)}
    write_json(os.path.join(d, "kappa.json"), kappa)
    log(f"adjudicate: {name}: kappa {kappa}")
    if todo and not finalize:
        write_csv(adj_path, todo)
        log(f"adjudicate: {len(todo)} disagreements written to {adj_path}; after the third expert fills the "
            f"final_* columns, rerun with --finalize")
        return adj_path
    path = os.path.join(d, "final_labels.csv")
    write_csv(path, out)
    log(f"adjudicate: final labels for {len(out)} sampled entries written")
    return path


def step_audit_sample(cfg: dict, cohort_file: str, labels_file: str, n_each: int = 50) -> str:
    """X3: random adjudicated negatives, 50 intermetallic and 50 other, for a third expert."""
    name, d = _dirs(cfg, labels_file)
    rows = {r["pred_id"]: r for r in read_csv(cohort_file)}
    fin = read_csv(os.path.join(d, "final_labels.csv"))
    neg = [f for f in fin if not f["final_first_report_date"]]
    inter = [f for f in neg if rows[f["pred_id"]]["chemistry_class"] == "intermetallic"]
    other = [f for f in neg if rows[f["pred_id"]]["chemistry_class"] != "intermetallic"]
    rng = np.random.default_rng(int(cfg["seeds"]["residual_audit"]))
    pick = lambda pool: [pool[i] for i in rng.choice(len(pool), size=min(n_each, len(pool)), replace=False)] if pool else []
    chosen = [("intermetallic", x) for x in pick(inter)] + [("other", x) for x in pick(other)]
    sheet = [{"audit_id": f"R{j:04d}", "stratum": s, "formula": rows[x["pred_id"]]["formula"],
              "elements": rows[x["pred_id"]]["elements"], "entry_date": rows[x["pred_id"]]["entry_date"],
              "found_missed_report": "", "report_date": "", "reference": "", "notes": ""}
             for j, (s, x) in enumerate(chosen)]
    keyrows = [{"audit_id": s["audit_id"], "pred_id": x["pred_id"], "stratum": st}
               for s, (st, x) in zip(sheet, chosen)]
    write_csv(os.path.join(d, "residual_audit_sheet.csv"), sheet)
    write_csv(os.path.join(d, "residual_audit_key.csv"), keyrows)
    write_json(os.path.join(d, "residual_audit_pool.json"),
               {"negatives_intermetallic": len(inter), "negatives_other": len(other)})
    log(f"audit-sample: {len(sheet)} adjudicated negatives for the third expert ({len(inter)} intermetallic, "
        f"{len(other)} other in the pool)")
    return os.path.join(d, "residual_audit_sheet.csv")
