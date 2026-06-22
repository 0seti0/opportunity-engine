"""WEEKLY RUN — runs the engine on your Claude Max SUBSCRIPTION (no API key). Cron this.

Flow:  discovery (deterministic, upstream) -> evidence bundle (deterministic) -> claude -p judgment
       (subscription) -> fact-check OVERRIDE (deterministic) -> rank -> snapshot + week-over-week delta.

The discovery/.py stages (attribute, cysdb, structural_cys, momentum, factcheck) refresh the candidate
list + artifacts upstream; this orchestrates the per-target judgment and the weekly digest.

Schedule (Mondays 8am) — `crontab -e`:
  0 8 * * 1 cd /path/to/layer1 && ~/.local/bin/uv run --python 3.12 python run_weekly.py >> runs/cron.log 2>&1
"""
import datetime
import json
from pathlib import Path

from assemble_evidence import assemble_evidence
from llm_claude_code import judge

HERE = Path(__file__).parent
RUNS = HERE / "runs"
RUNS.mkdir(exist_ok=True)

VERDICT_SCHEMA = {
    "type": "object",
    "properties": {
        "target": {"type": "string"},
        "verdict": {"type": "string", "enum": ["STRONG_BIC", "POSSIBLE_BIC", "WEAK", "DEAD"]},
        "lane": {"type": "string", "enum": ["OPEN", "CONTESTED", "CLOSED"]},
        "score": {"type": "number", "description": "0-100 overall covalent-BIC attractiveness"},
        "thesis": {"type": "string", "description": "one-line covalent best-in-class thesis"},
        "top_risk": {"type": "string"},
    },
    "required": ["target", "verdict", "lane", "score", "thesis", "top_risk"],
}


def candidates():
    # discovery output -> targets to judge this week. (Swap for the live discovery list; here: shortlist.)
    u = json.load(open(HERE / "unified_candidates.json"))
    return [r["sym"] for r in u if r.get("verdict") in ("HOT_OPEN", "OPEN", "POSSIBLE")]


def _delta(prev_path, cur_rows):
    prev = {r["target"]: r for r in json.load(open(prev_path))}
    cur = {r["target"]: r for r in cur_rows}
    out = []
    for t, r in cur.items():
        if t not in prev:
            out.append(f"NEW  {t} ({r['verdict']}, score {r.get('score')})")
        elif r.get("verdict") != prev[t].get("verdict"):
            out.append(f"Δ    {t}: {prev[t].get('verdict')} -> {r.get('verdict')}")
    out += [f"gone {t}" for t in prev if t not in cur]
    return out or ["no changes vs last week"]


def run():
    rows = []
    for sym in candidates():
        ev = assemble_evidence(sym)
        try:
            v = judge(ev, VERDICT_SCHEMA)
        except Exception as e:
            v = {"target": sym, "verdict": "ERROR", "lane": "?", "score": 0,
                 "thesis": f"judgment failed: {str(e)[:120]}", "top_risk": ""}
        if ev.get("facts"):                       # deterministic fact-check OVERRIDES the model
            v["selectivity"] = ev["facts"].get("selectivity")
        rows.append(v)
    rows.sort(key=lambda r: -r.get("score", 0))

    today = datetime.date.today().isoformat()     # plain script (not a workflow) -> datetime is fine
    json.dump(rows, open(RUNS / f"run_{today}.json", "w"), indent=1)
    prev = sorted(RUNS.glob("run_*.json"))
    delta = _delta(prev[-2], rows) if len(prev) >= 2 else ["first run — no baseline"]

    print(f"=== weekly run {today}: {len(rows)} targets (judged on Max subscription, no API key) ===")
    for i, r in enumerate(rows[:15], 1):
        print(f"{i:2}. {r['target']:9} {r['verdict']:12} score={int(r.get('score', 0)):>3} "
              f"sel={str(r.get('selectivity', '?')):<20} {r['lane']:9} {r['thesis'][:55]}")
    print("\nΔ vs last week:")
    for d in delta:
        print("   " + d)


if __name__ == "__main__":
    run()
