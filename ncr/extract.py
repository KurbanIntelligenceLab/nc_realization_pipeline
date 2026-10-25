"""Language-model extraction, the grounding check, entry-level model labels and freezing.

The model extracts facts from each retrieved document (Appendix A and B of the
runbook); this code, not the model, assigns tiers and dates. A positive record
counts only if its quoted passage occurs in the retrieved text (normalized, at
least 40 characters); otherwise it is an abstention (Methods).
"""
from __future__ import annotations

import hashlib
import json
import os
import re
import sys
import time

from . import config as C
from .util import date_str, log, parse_date, read_csv, read_json, sha256_file, write_csv, write_json

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)
import nc_realization_cohorts as NCC                 # noqa: E402
import nc_realization_metrics as NCM                 # noqa: E402

SYSTEM = ("You extract experimental synthesis reports for inorganic compounds. Answer only from the "
          "supplied document. If the document does not clearly report an experimental synthesis of "
          "the target, answer \"no\" or \"abstain\". Never guess a year.")

USER_TEMPLATE = ("Target formula: {formula}\nElement set: {elements}\n"
                 "Predicted structure (if structure-level): space group {sg}, prototype {proto}\n"
                 "Document metadata: {doi}, {publication_date}\nDocument text: {text}\n"
                 "Return JSON matching the schema.")

_YNU = {"type": "string", "enum": ["yes", "no", "unclear"]}
_NSTR = {"type": ["string", "null"]}
SCHEMA = {
    "type": "object", "additionalProperties": False,
    "properties": {
        "doi": {"type": "string"},
        "reports_synthesis": {"type": "string", "enum": ["yes", "no", "abstain"]},
        "formula_as_reported": {"type": "string"},
        "composition_range": {"type": "object", "additionalProperties": False,
                              "properties": {"end_a": _NSTR, "end_b": _NSTR}, "required": ["end_a", "end_b"]},
        "lattice_parameters": _NSTR,
        "composition_matches_target": _YNU,
        "structure_reported": {"type": "object", "additionalProperties": False,
                               "properties": {"space_group": _NSTR, "prototype": _NSTR, "partial_occupancy": _YNU},
                               "required": ["space_group", "prototype", "partial_occupancy"]},
        "conditions": {"type": "object", "additionalProperties": False,
                       "properties": {"high_pressure": _YNU, "thin_film_only": _YNU},
                       "required": ["high_pressure", "thin_film_only"]},
        "publication_date": {"type": "string"},
        "received_date": _NSTR,
        "evidence_passage": {"type": "string"},
        "evidence_location": {"type": "string"},
        "confidence": {"type": "string", "enum": ["high", "medium", "low"]},
    },
    "required": ["doi", "reports_synthesis", "formula_as_reported", "composition_range", "lattice_parameters",
                 "composition_matches_target", "structure_reported", "conditions", "publication_date",
                 "received_date", "evidence_passage", "evidence_location", "confidence"],
}


def prompt_hash(model: str) -> str:
    blob = json.dumps({"system": SYSTEM, "user": USER_TEMPLATE, "schema": SCHEMA, "model": model}, sort_keys=True)
    return hashlib.sha256(blob.encode()).hexdigest()


# ----------------------------------------------------------------------------
# Providers
# ----------------------------------------------------------------------------
class AnthropicProvider:
    def __init__(self, model: str, max_tokens: int = 1024, temperature=None):
        try:
            import anthropic
        except ImportError as exc:
            raise ImportError("pip install anthropic (see requirements-pipeline.txt)") from exc
        if not os.environ.get("ANTHROPIC_API_KEY"):
            raise RuntimeError("set ANTHROPIC_API_KEY in your terminal before running 'extract'")
        self.client = anthropic.Anthropic()
        self.model, self.max_tokens, self.temperature = model, max_tokens, temperature
        self.tokens_in = self.tokens_out = 0

    def __call__(self, system: str, user: str) -> dict:
        kw = {}
        if self.temperature is not None:
            kw["extra_body"] = {"temperature": self.temperature}
        for attempt in range(6):
            try:
                msg = self.client.messages.create(
                    model=self.model, max_tokens=self.max_tokens, system=system,
                    messages=[{"role": "user", "content": user}],
                    output_config={"format": {"type": "json_schema", "schema": SCHEMA}}, **kw)
                self.tokens_in += msg.usage.input_tokens
                self.tokens_out += msg.usage.output_tokens
                text = "".join(getattr(b, "text", "") for b in msg.content)
                return json.loads(text)
            except json.JSONDecodeError:
                return {"_error": "unparseable model output"}
            except Exception as exc:  # noqa: BLE001
                wait = min(60, 2 ** attempt)
                log(f"extract: API error {type(exc).__name__}: {exc}; retrying in {wait} s")
                time.sleep(wait)
        return {"_error": "API failed repeatedly"}


class MockProvider:
    """Deterministic stand-in for tests and dry runs. It answers 'yes' when a sentence
    contains the target formula and a synthesis verb, and injects reproducible errors
    (missed reports and ungrounded claims) so the downstream checks are exercised."""
    VERBS = ("synthesized", "synthesised", "prepared", "grown", "obtained")

    def __init__(self, miss_rate: float = 0.1, fake_rate: float = 0.03):
        self.miss_rate, self.fake_rate = miss_rate, fake_rate
        self.tokens_in = self.tokens_out = 0

    @staticmethod
    def _u(s: str) -> float:
        return int(hashlib.sha256(s.encode()).hexdigest()[:8], 16) / 0xFFFFFFFF

    def __call__(self, system: str, user: str) -> dict:
        f_in = re.search(r"Target formula: (\S+)", user).group(1)
        from pymatgen.core import Composition
        f = Composition(f_in).reduced_formula
        meta = re.search(r"Document metadata: (.*?), (\S*)\n", user)
        text = user.split("Document text: ", 1)[1].rsplit("\nReturn JSON", 1)[0]
        base = {"doi": meta.group(1) if meta else "", "formula_as_reported": f,
                "composition_range": {"end_a": None, "end_b": None}, "lattice_parameters": None,
                "composition_matches_target": "yes",
                "structure_reported": {"space_group": None, "prototype": None, "partial_occupancy": "unclear"},
                "conditions": {"high_pressure": "no", "thin_film_only": "no"},
                "publication_date": meta.group(2) if meta else "", "received_date": None,
                "evidence_location": "abstract", "confidence": "high"}
        sg = re.search(r"space group (\d+)", text)
        if sg:
            base["structure_reported"] = {"space_group": sg.group(1), "prototype": None, "partial_occupancy": "no"}
        u = self._u(user)
        for sent in re.split(r"(?<=\.)\s+", text):
            if f in sent and any(v in sent for v in self.VERBS):
                if u < self.miss_rate:
                    break
                return dict(base, reports_synthesis="yes", evidence_passage=sent.strip())
        if u > 1 - self.fake_rate:
            return dict(base, reports_synthesis="yes",
                        evidence_passage=f"We report the first synthesis of {f} by a solid-state route at 900 C.")
        return dict(base, reports_synthesis="no", evidence_passage="")


def make_provider(cfg: dict, kind: str):
    if kind == "mock":
        return MockProvider()
    lc = cfg.get("llm", {})
    return AnthropicProvider(C.need(cfg, "llm.model"), int(lc.get("max_tokens", 1024)), lc.get("temperature"))


# ----------------------------------------------------------------------------
# Record-level rules (applied identically to expert records where relevant)
# ----------------------------------------------------------------------------
def _sg_number(sym):
    if sym in (None, ""):
        return None
    try:
        return int(sym)
    except (TypeError, ValueError):
        pass
    try:
        from pymatgen.symmetry.groups import SpaceGroup
        return SpaceGroup(str(sym).replace(" ", "")).int_number
    except Exception:  # noqa: BLE001
        return None


def record_tier(rec: dict, target: dict, comp_delta: float) -> str | None:
    """Tier of a grounded positive literature record. Composition rule as in
    Table 2. Without a refined structure, structure agreement is judged from the
    reported space group or prototype and flagged (a model-side proxy; experts
    decide T1/T2 for the certified labels)."""
    rep = rec.get("formula_as_reported") or ""
    comp_ok = False
    try:
        comp_ok = NCC.composition_distance(rep, target["formula"]) <= comp_delta
    except Exception:  # noqa: BLE001
        comp_ok = False
    rng = rec.get("composition_range") or {}
    if not comp_ok and rng.get("end_a") and rng.get("end_b"):
        try:
            comp_ok = NCC.in_solid_solution_range(target["formula"], rng["end_a"], rng["end_b"], comp_delta)
        except Exception:  # noqa: BLE001
            comp_ok = False
    if not comp_ok:
        return None
    if target.get("granularity") == "structure":
        st = rec.get("structure_reported") or {}
        sg_t = _sg_number(target.get("spacegroup_number"))
        sg_r = _sg_number(st.get("space_group"))
        proto_t = (target.get("prototype") or "").strip().lower()
        proto_r = (st.get("prototype") or "").strip().lower()
        if (sg_t and sg_r and sg_t == sg_r) or (proto_t and proto_r and proto_t == proto_r):
            return "T2"
    return "T3"


def _primary(gran: str) -> str:
    return "T2" if gran == "structure" else "T3"


# ----------------------------------------------------------------------------
# Step: extract
# ----------------------------------------------------------------------------
def step_extract(cfg: dict, cohort_file: str, provider_kind: str, cache_name: str = "extracted") -> str:
    rows = read_csv(cohort_file)
    prov = make_provider(cfg, provider_kind)
    model = "mock" if provider_kind == "mock" else C.need(cfg, "llm.model")
    ph = prompt_hash(model)
    rdir = os.path.join(cfg["workdir"], "retrieved")
    edir = os.path.join(cfg["workdir"], cache_name)
    os.makedirs(edir, exist_ok=True)
    n_calls = 0
    for i, r in enumerate(rows):
        from .retrieve import _safe
        out = os.path.join(edir, _safe(r["pred_id"]) + ".json")
        if os.path.exists(out):
            continue
        rp = os.path.join(rdir, _safe(r["pred_id"]) + ".json")
        if not os.path.exists(rp):
            raise FileNotFoundError(f"run 'retrieve' first ({rp} missing)")
        docs = read_json(rp)["docs"]
        recs = []
        for d in docs:
            user = USER_TEMPLATE.format(formula=r["formula"], elements=r["elements"],
                                        sg=r.get("spacegroup_number") or "n/a", proto=r.get("prototype") or "n/a",
                                        doi=d.get("doi") or d["doc_id"], publication_date=d.get("publication_date") or "",
                                        text=d["text"])
            rec = prov(SYSTEM, user)
            n_calls += 1
            rec["_doc_id"] = d["doc_id"]
            rec["_doc_publication_date"] = d.get("publication_date") or ""
            rec["_doc_received_date"] = d.get("received_date") or ""
            rec["_grounded"] = bool(rec.get("reports_synthesis") == "yes"
                                    and NCM.passage_in_document(rec.get("evidence_passage", ""), d["text"]))
            recs.append(rec)
        write_json(out, {"pred_id": r["pred_id"], "model": model, "prompt_hash": ph, "records": recs})
        if (i + 1) % 100 == 0:
            log(f"extract: {i + 1}/{len(rows)} entries, {n_calls} calls, tokens in/out {prov.tokens_in}/{prov.tokens_out}")
    log(f"extract: {os.path.basename(cohort_file)} done; calls {n_calls}; tokens in/out {prov.tokens_in}/{prov.tokens_out}")
    cost = os.path.join(cfg["workdir"], "cost_log.csv")
    prev = read_csv(cost) if os.path.exists(cost) else []
    prev.append({"cohort_file": os.path.basename(cohort_file), "cache": cache_name, "model": model,
                 "prompt_hash": ph, "calls": n_calls, "tokens_in": prov.tokens_in, "tokens_out": prov.tokens_out,
                 "time": time.strftime("%Y-%m-%d %H:%M:%S")})
    write_csv(cost, prev)
    return edir


# ----------------------------------------------------------------------------
# Step: label (entry-level model labels from records and database matches)
# ----------------------------------------------------------------------------
def entry_label(r: dict, recs: list[dict], db: dict | None, comp_delta: float) -> dict:
    gran = r.get("granularity", "structure")
    prim = _primary(gran)
    dates, rdates, flags, abstain = [], [], set(), False
    for rec in recs:
        if rec.get("_error"):
            abstain = True
            continue
        if rec.get("reports_synthesis") == "abstain":
            abstain = True
            continue
        if rec.get("reports_synthesis") != "yes":
            continue
        if not rec.get("_grounded"):
            abstain = True                                     # ungrounded positive = abstention
            flags.add("ungrounded_claim")
            continue
        tier = record_tier(rec, r, comp_delta)
        if tier is None:
            if rec.get("composition_matches_target") == "unclear":
                abstain = True
            continue
        if not NCC.meets(tier, prim):
            continue
        d, _ = parse_date(rec.get("_doc_publication_date") or rec.get("publication_date"))
        if d is None:
            abstain = True
            flags.add("undated_report")
            continue
        dates.append(d)
        rd, _ = parse_date(rec.get("received_date") or rec.get("_doc_received_date"))
        rdates.append(rd or d)
        if gran == "structure":
            flags.add("structure_from_metadata")
        cond = rec.get("conditions") or {}
        if cond.get("high_pressure") == "yes":
            flags.add("high_pressure")
        if cond.get("thin_film_only") == "yes":
            flags.add("thin_film_only")
    if db:
        d, _ = parse_date(db.get(f"db_first_{prim}", ""))
        if d is not None:
            dates.append(d)
            rdates.append(d)
            flags.update(x for x in (db.get("db_flags") or "").split() if x)
    if dates:
        label, first, first_r = "report", min(dates), min(rdates)
    else:
        label, first, first_r = ("abstain" if abstain else "none"), None, None
    return {"model_label": label, "model_first_date": date_str(first),
            "model_first_received_date": date_str(first_r), "model_flags": " ".join(sorted(flags))}


def first_mention(formula: str, docs: list[dict]) -> str:
    """Earliest retrieved document whose text contains the reduced formula
    (attention outcome, X9). Model-side proxy; experts record the certified value."""
    from pymatgen.core import Composition
    red = Composition(formula).reduced_formula
    pat = re.compile(r"(?<![A-Za-z0-9])" + re.escape(red) + r"(?![a-z0-9])")
    ds = [parse_date(d.get("publication_date"))[0] for d in docs if pat.search(d.get("text", ""))]
    ds = [d for d in ds if d]
    return date_str(min(ds)) if ds else ""


def step_label(cfg: dict, cohort_file: str, cache_name: str = "extracted", out_path: str | None = None,
               db_path: str | None = None) -> str:
    from .retrieve import _safe
    rows = read_csv(cohort_file)
    comp_delta = float((cfg.get("identity_tolerances") or {}).get("comp_delta", 0.02))
    dbp = db_path or cohort_file.replace(os.sep + "cohorts" + os.sep, os.sep + "dbmatch" + os.sep)
    db = {x["pred_id"]: x for x in read_csv(dbp)} if os.path.exists(dbp) else {}
    if not db:
        log(f"label: no database matches for {os.path.basename(cohort_file)} (run 'dbmatch' to include them)")
    out = []
    for r in rows:
        e = read_json(os.path.join(cfg["workdir"], cache_name, _safe(r["pred_id"]) + ".json"))
        docs = read_json(os.path.join(cfg["workdir"], "retrieved", _safe(r["pred_id"]) + ".json"))["docs"]
        lab = entry_label(r, e["records"], db.get(r["pred_id"]), comp_delta)
        out.append({"pred_id": r["pred_id"], "cohort": r["cohort"], **lab,
                    "model_first_mention_date": first_mention(r["formula"], docs),
                    "prompt_hash": e["prompt_hash"], "model": e["model"]})
    sub = "labels" if cache_name == "extracted" else f"labels_{cache_name}"
    path = out_path or os.path.join(cfg["workdir"], sub, os.path.basename(cohort_file))
    if os.path.exists(path + ".sha256"):
        raise RuntimeError(f"{path} is frozen; labels cannot be regenerated after freezing")
    write_csv(path, out)
    cnt = {k: sum(1 for o in out if o["model_label"] == k) for k in ("report", "none", "abstain")}
    log(f"label: {os.path.basename(cohort_file)}: {cnt}")
    return path


def step_freeze(cfg: dict, labels_file: str, manifest) -> str:
    h = sha256_file(labels_file)
    p = labels_file + ".sha256"
    if os.path.exists(p):
        old = open(p, encoding="utf-8").read().strip()
        if old != h:
            raise RuntimeError(f"{labels_file} changed after it was frozen (frozen {old[:12]}, now {h[:12]}); "
                               "labels must not change after freezing")
        log(f"freeze: already frozen ({h[:12]})")
        return h
    with open(p, "w", encoding="utf-8") as fh:
        fh.write(h + "\n")
    os.chmod(labels_file, 0o444)
    manifest.add("freeze", labels_file=os.path.relpath(labels_file, cfg["workdir"]), sha256=h)
    log(f"freeze: {os.path.basename(labels_file)} sha256 {h}")
    return h
