"""Literature retrieval for each prediction.

Two backends:
  local     a JSONL corpus you provide (one document per line: doc_id, doi, title,
            text, publication_date, received_date). Use this for licensed full text.
  openalex  the public OpenAlex works search (titles and abstracts). Requires a
            contact e-mail (OpenAlex "polite pool"); results are cached per entry.
Queries follow Table S1: reduced formula, element-order permutations, a fractional
notation, and the element set.
"""
from __future__ import annotations

import itertools
import json
import os
import re
import time
import urllib.parse
import urllib.request

from .util import log, parse_date, read_csv, write_json

MAX_PERMUTATION_ELEMENTS = 4


def query_variants(formula: str) -> list[str]:
    from pymatgen.core import Composition
    c = Composition(formula)
    red = c.reduced_formula
    amounts = c.reduced_composition.get_el_amt_dict()
    els = list(amounts)
    out = [red]
    if len(els) <= MAX_PERMUTATION_ELEMENTS:
        for perm in itertools.permutations(els):
            s = "".join(e + (_num(amounts[e])) for e in perm)
            if s not in out:
                out.append(s)
    frac = c.fractional_composition.get_el_amt_dict()
    out.append("".join(f"{e}{frac[e]:.3g}" for e in els))
    out.append(" ".join(sorted(els)))                     # element set
    return out


def _num(x: float) -> str:
    return "" if abs(x - 1) < 1e-9 else (str(int(round(x))) if abs(x - round(x)) < 1e-9 else f"{x:g}")


# ----------------------------------------------------------------------------
class LocalCorpus:
    def __init__(self, path: str):
        self.docs = []
        with open(path, encoding="utf-8") as fh:
            for ln in fh:
                if ln.strip():
                    d = json.loads(ln)
                    for k in ("doc_id", "text", "publication_date"):
                        if k not in d:
                            raise ValueError(f"corpus document missing '{k}'")
                    self.docs.append(d)

    def search(self, variants: list[str], k: int, until: str) -> list[dict]:
        until_d, _ = parse_date(until)
        pats = [re.compile(r"(?<![A-Za-z0-9])" + re.escape(v) + r"(?![a-z0-9])") for v in variants[:-1]]
        hits = []
        for d in self.docs:
            pd_, _ = parse_date(d["publication_date"])
            if pd_ and until_d and pd_ > until_d:
                continue
            txt = (d.get("title") or "") + "\n" + d["text"]
            if any(p.search(txt) for p in pats):
                hits.append({"doc_id": d["doc_id"], "doi": d.get("doi", ""), "title": d.get("title", ""),
                             "text": txt, "publication_date": d["publication_date"],
                             "received_date": d.get("received_date", ""), "backend": "local"})
        hits.sort(key=lambda x: x["publication_date"])
        return hits[:k]


class OpenAlex:
    URL = "https://api.openalex.org/works"

    def __init__(self, email: str, pause: float = 0.12):
        if not email or "@" not in email:
            raise ValueError("OpenAlex needs a contact e-mail (config retrieval.openalex_email)")
        self.email, self.pause = email, pause

    def _get(self, params: dict) -> dict:
        params = dict(params, mailto=self.email)
        url = self.URL + "?" + urllib.parse.urlencode(params)
        for attempt in range(5):
            try:
                with urllib.request.urlopen(url, timeout=60) as fh:
                    time.sleep(self.pause)
                    return json.loads(fh.read().decode("utf-8"))
            except Exception as exc:  # noqa: BLE001
                wait = 2 ** attempt
                log(f"openalex: {exc}; retrying in {wait} s")
                time.sleep(wait)
        raise RuntimeError(f"OpenAlex request failed repeatedly: {url}")

    @staticmethod
    def _abstract(inv: dict | None) -> str:
        if not inv:
            return ""
        pos = []
        for w, idx in inv.items():
            pos.extend((i, w) for i in idx)
        return " ".join(w for _, w in sorted(pos))

    def search(self, variants: list[str], k: int, until: str) -> list[dict]:
        seen, hits = set(), []
        for v in variants:
            data = self._get({"search": f"\"{v}\"", "per-page": min(k, 50),
                              "filter": f"to_publication_date:{until}",
                              "select": "id,doi,title,publication_date,abstract_inverted_index"})
            for w in data.get("results", []):
                if w["id"] in seen:
                    continue
                seen.add(w["id"])
                txt = (w.get("title") or "") + "\n" + self._abstract(w.get("abstract_inverted_index"))
                hits.append({"doc_id": w["id"], "doi": w.get("doi") or "", "title": w.get("title") or "",
                             "text": txt, "publication_date": w.get("publication_date") or "",
                             "received_date": "", "backend": "openalex"})
            if len(hits) >= k:
                break
        hits.sort(key=lambda x: x["publication_date"] or "9999")
        return hits[:k]


def make_backend(cfg: dict, path_resolver):
    rc = cfg.get("retrieval", {})
    kind = rc.get("backend", "local")
    if kind == "local":
        return LocalCorpus(path_resolver(rc["local_corpus"]))
    if kind == "openalex":
        return OpenAlex(rc.get("openalex_email", ""))
    raise ValueError(f"unknown retrieval backend {kind!r}")


def step_retrieve(cfg: dict, cohort_file: str, path_resolver, cutoff_date: str) -> str:
    rows = read_csv(cohort_file)
    k = int(cfg.get("retrieval", {}).get("k", 20))
    backend = make_backend(cfg, path_resolver)
    outdir = os.path.join(cfg["workdir"], "retrieved")
    os.makedirs(outdir, exist_ok=True)
    n_new = 0
    for i, r in enumerate(rows):
        p = os.path.join(outdir, _safe(r["pred_id"]) + ".json")
        if os.path.exists(p):
            continue
        docs = backend.search(query_variants(r["formula"]), k, cutoff_date)
        write_json(p, {"pred_id": r["pred_id"], "formula": r["formula"], "docs": docs})
        n_new += 1
        if (i + 1) % 200 == 0:
            log(f"retrieve: {i + 1}/{len(rows)}")
    log(f"retrieve: {os.path.basename(cohort_file)}: {n_new} new entries searched ({len(rows)} total)")
    return outdir


def _safe(pred_id: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", pred_id)
