"""Emerging-target query: rank targets from the claim log by earliest-signal recency + low crowding.
Corroboration (count of evidence types) is a confidence TAG, NOT the rank driver -- the whole point is to
surface niche/early targets *before* they're validated and crowded, which corroboration-as-rank would bury."""
import json
from collections import defaultdict
from pathlib import Path

LOG = Path(__file__).parent / "claims.jsonl"


def rank():
    by_target = defaultdict(list)
    for line in LOG.read_text().splitlines():
        c = json.loads(line)
        by_target[c["subject"]].append(c)
    rows = []
    for hid, cs in by_target.items():
        actors = {c["actor"] for c in cs if c["actor"]}
        dates = [c["valid_time"] for c in cs if c["valid_time"]]
        if not actors or not dates:              # quality floor: a credible, dated signal. ponytail: Path A precision is the real noise filter, upstream
            continue
        rows.append({"target": hid, "earliest": min(dates), "sponsors": len(actors),
                     "evidence": sorted({c["evidence_type"] for c in cs}), "signals": len(cs)})
    rows.sort(key=lambda r: (r["earliest"], -r["sponsors"]), reverse=True)   # emerging-first: recent first-signal, then least crowded
    return rows


if __name__ == "__main__":
    r = rank()
    assert all(r[i]["earliest"] >= r[i + 1]["earliest"] for i in range(len(r) - 1))   # emerging-first ordering holds
    for x in r[:12]:
        tag = "+".join(x["evidence"]) + (" ✓corroborated" if len(x["evidence"]) > 1 else "")
        print(f"  {x['target']:11} earliest {x['earliest']}  {x['sponsors']:2} sponsor(s)  [{tag}]")
    print(f"{len(r)} targets ranked emerging-first")
