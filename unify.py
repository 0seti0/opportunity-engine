"""Step 3 — UNIFY Leg 1 (patent) + Leg 2 (clinical) into one covalent-move log, keyed by company + target,
dated by the EARLIEST signal. The two legs have opposite timing: a CLINICAL signal is early but carries no
chemistry; a PATENT is late but carries the warhead. Cross-linking them by company+target gives the holy
grail — an early-dated covalent move WITH the chemistry once the patent lands. -> the alert digest (Step 4).

  python unify.py            # build + print the covalent-move digest (cross-linked, ranked, diffed vs last run)
  python unify.py --selftest # offline check of the merge/cross-link logic
"""
import datetime
import json
import sys
from collections import defaultdict
from pathlib import Path

from attribute import trusted
from whos_moving import firm

HERE = Path(__file__).parent
SEEN = HERE / ".unify_seen.json"
COVALENT = ("LIKELY", "FINGERPRINT", "POSSIBLE")     # clinical covalency tiers worth logging


def _patent_moves():
    tg = json.load(open(HERE / "patent_targets.json"))
    feed = {e[0]: e for e in json.load(open(HERE / "covalent_warhead_feed_clean.json"))}
    for pn, d in tg.items():
        e = feed.get(pn)
        if d.get("sym") and trusted(d) and e:
            yield {"company": firm(e[2] or ""), "target": d["sym"], "date": e[1] or "", "type": "patent",
                   "id": pn, "chemistry": e[5] if len(e) > 5 else [], "stage": "patent"}


def _clinical_moves():
    p = HERE / "clinical_events.json"
    for ev in (json.load(open(p)) if p.exists() else []):
        if ev.get("covalency") in COVALENT:
            yield {"company": firm(ev.get("sponsor", "")), "target": ev.get("target", ""),
                   "date": ev.get("posted", ""), "type": "clinical", "id": ev.get("nct", ""),
                   "chemistry": [], "stage": ev.get("covalency")}


def _unify(moves):
    """Combined [move,...] -> one entry per (company, target), earliest-dated, cross-link flagged."""
    by = defaultdict(lambda: {"signals": [], "chemistry": []})
    for m in moves:
        if m.get("company") and m.get("target"):
            by[(m["company"], m["target"])]["signals"].append(m)
            by[(m["company"], m["target"])]["chemistry"] += m.get("chemistry", [])
    out = []
    for (company, target), v in by.items():
        sigs = sorted(v["signals"], key=lambda s: s["date"] or "9999")
        types = {s["type"] for s in sigs}
        out.append({"company": company, "target": target,
                    "earliest": sigs[0]["date"], "earliest_type": sigs[0]["type"],
                    "cross_linked": len(types) > 1, "has_chemistry": bool(v["chemistry"]),
                    "warheads": sorted(set(v["chemistry"])),
                    "n_patent": sum(s["type"] == "patent" for s in sigs),
                    "n_clinical": sum(s["type"] == "clinical" for s in sigs),
                    "signals": [{"type": s["type"], "id": s["id"], "date": s["date"]} for s in sigs]})
    # rank: cross-linked first (both signals), then clinical-first (earliest warning), then has-chemistry, then recency
    out.sort(key=lambda e: (not e["cross_linked"], e["earliest_type"] != "clinical",
                            not e["has_chemistry"], e["earliest"] < "0", e["earliest"]), reverse=False)
    out.sort(key=lambda e: e["earliest"], reverse=True)            # recency within tier
    out.sort(key=lambda e: (not e["cross_linked"], e["earliest_type"] != "clinical"))
    return out


def unify():
    return _unify(list(_patent_moves()) + list(_clinical_moves()))


def run():
    log = unify()
    today = datetime.date.today().isoformat()
    xlinked = [e for e in log if e["cross_linked"]]
    clinical = [e for e in log if e["n_clinical"] and not e["cross_linked"]]
    seen = {tuple(x) for x in json.loads(SEEN.read_text())} if SEEN.exists() else set()
    new = [e for e in log if (e["company"], e["target"]) not in seen]

    print(f"=== covalent-move log {today}: {len(log)} company×target moves "
          f"({len(xlinked)} cross-linked · {len(clinical)} clinical-only) ===")
    print(f"\n⚡ CROSS-LINKED — clinical + patent on the same company×target (early signal WITH chemistry): {len(xlinked)}")
    for e in xlinked[:15]:
        print(f"   {e['company'][:22]:22} {e['target']:8} earliest {e['earliest']} ({e['earliest_type']}) "
              f"· {e['n_patent']}pat/{e['n_clinical']}clin · warheads={e['warheads']}")
    print(f"\n🧪 CLINICAL-ONLY — early warning, no patent yet (watch for the filing): {len(clinical)}")
    for e in clinical[:15]:
        ncts = [s["id"] for s in e["signals"] if s["type"] == "clinical"][:2]
        print(f"   {e['company'][:22]:22} {e['target']:8} {e['earliest']} {ncts}")
    print(f"\n🆕 NEW vs last run: {len(new)}")
    for e in new[:20]:
        print(f"   {e['company'][:22]:22} {e['target']:8} {e['earliest']} ({e['earliest_type']}) "
              f"{'⚗ ' + ','.join(e['warheads']) if e['has_chemistry'] else ''}")

    SEEN.write_text(json.dumps(sorted([[e["company"], e["target"]] for e in log])))
    json.dump(log, open(HERE / "covalent_move_log.json", "w"), indent=1)
    print(f"\nsaved covalent_move_log.json ({len(log)} moves)")
    return log


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        moves = [{"company": "ACME", "target": "EGFR", "type": "patent", "date": "2025-03-01",
                  "id": "WO-1", "chemistry": ["acrylamide"]},
                 {"company": "ACME", "target": "EGFR", "type": "clinical", "date": "2024-06-01",
                  "id": "NCT-1", "chemistry": []},
                 {"company": "BETA", "target": "BTK", "type": "patent", "date": "2026-01-01",
                  "id": "WO-2", "chemistry": ["butynamide"]}]
        log = _unify(moves)
        acme = next(e for e in log if e["company"] == "ACME")
        assert acme["cross_linked"] and acme["has_chemistry"], acme           # both legs on ACME×EGFR
        assert acme["earliest"] == "2024-06-01" and acme["earliest_type"] == "clinical"  # earliest-dated by the trial
        assert log[0]["company"] == "ACME"                                    # cross-linked ranks first
        beta = next(e for e in log if e["company"] == "BETA")
        assert not beta["cross_linked"] and beta["warheads"] == ["butynamide"]
        print("unify self-check OK:", [(e["company"], e["target"], e["cross_linked"]) for e in log])
    else:
        run()
