"""Deterministic evidence bundle per target — gathered from the engine's existing artifacts + live
ClinicalTrials.gov. NO LLM. This is the structured evidence the judgment reasons over; `facts` (from
the deterministic UniProt fact-check) is the ground truth that OVERRIDES the model afterward.
"""
import json
from pathlib import Path

import backtest as bt   # ct_trials (stdlib urllib + cache)

HERE = Path(__file__).parent


def _load(name):
    p = HERE / name
    return json.load(open(p)) if p.exists() else None


def assemble_evidence(sym):
    ev = {"target": sym}
    # covalent handle: CysDB ligandable cysteines + KLIFS pocket-cysteine class
    cy = {r["sym"]: r for r in (_load("cysdb_opportunities.json") or {}).get("opportunities", [])}.get(sym)
    st = {r["sym"]: r for r in (_load("structural_cys_opportunities.json") or {}).get("opportunities", [])}.get(sym)
    ev["cysdb_ligandable_cys"] = cy.get("n_lig_cys") if cy else None
    ev["pocket_cys_handles"] = st.get("handles") if st else None
    # GROUND-TRUTH facts (residue identity + paralog selectivity) — overrides the model later
    ev["facts"] = (_load("factcheck.json") or {}).get(sym)
    # validation: live ClinicalTrials.gov interventional counts
    try:
        ev["clinical_trials_total"] = bt.ct_trials(sym, "2000-01-01", "2026-12-31")
        ev["clinical_trials_recent_24_26"] = bt.ct_trials(sym, "2024-01-01", "2026-12-31")
    except Exception:
        ev["clinical_trials_total"] = ev["clinical_trials_recent_24_26"] = None
    # lane: covalent patents (title/claims only) attributed to this target
    tg = _load("patent_targets.json") or {}
    feed = {e[0]: e for e in (_load("covalent_warhead_feed_clean.json") or [])}
    pats = [{"patent": pn, "title": feed[pn][3], "date": feed[pn][1]}
            for pn, d in tg.items()
            if d.get("sym") == sym and d.get("conf") in ("title", "claims") and pn in feed]
    ev["covalent_patents"] = pats or "none in feed (lane likely open)"
    # momentum: literature acceleration (2025-26 vs 2021-23)
    scan = {r["sym"]: r for r in (_load("emerging_open_ranked.json") or {}).get("all", [])}.get(sym)
    ev["momentum_accel"] = scan.get("accel") if scan else None
    return ev


if __name__ == "__main__":
    print(json.dumps(assemble_evidence("SLC1A5"), indent=1))
