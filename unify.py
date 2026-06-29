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
from collections import Counter, defaultdict
from pathlib import Path

from attribute import trusted
from tractable import covalently_tractable
from whos_moving import firm

HERE = Path(__file__).parent
SEEN = HERE / ".unify_seen.json"
COVALENT = ("LIKELY", "FINGERPRINT")                 # confident clinical covalency only (POSSIBLE = noisy target-context)


def _patent_moves():
    tg = json.load(open(HERE / "patent_targets.json"))
    feed = {e[0]: e for e in json.load(open(HERE / "covalent_warhead_feed_clean.json"))}
    for pn, d in tg.items():
        e = feed.get(pn)
        if d.get("sym") and trusted(d) and covalently_tractable(d["sym"]) and e:    # drop antibody-class targets
            yield {"company": firm(e[2] or ""), "target": d["sym"], "date": e[1] or "", "type": "patent",
                   "id": pn, "chemistry": e[5] if len(e) > 5 else [], "stage": "patent"}


def _clinical_moves():
    p = HERE / "clinical_events.json"
    for ev in (json.load(open(p)) if p.exists() else []):
        if ev.get("covalency") in COVALENT and covalently_tractable(ev.get("target", "")):
            yield {"company": firm(ev.get("sponsor", "")), "target": ev.get("target", ""),
                   "date": ev.get("posted", ""), "type": "clinical", "id": ev.get("nct", ""),
                   "chemistry": [], "stage": ev.get("covalency")}


# mega-crowded covalent franchises — always noise for a fast-follower, regardless of our census's (sparse,
# Western-skewed) filer count. KRAS shows only 2 filers in our census but is the most saturated covalent lane.
SATURATED = {"KRAS", "EGFR", "BTK", "ERBB2", "KRASG12C"}


def _lane(target, n_filers):
    """Crowding tier — the Axiom-relevance axis. A move on a SATURATED franchise is noise; OPEN/EMERGING is the
    white space worth fast-following. Uses a curated denylist (reliable) OR the census filer count (data-driven)."""
    if target in SATURATED or n_filers > 10:
        return "CROWDED"
    return "OPEN" if n_filers <= 1 else "EMERGING" if n_filers <= 3 else "CONTESTED"


def _unify(moves):
    """Combined [move,...] -> one entry per (company, target), earliest-dated, cross-link + crowding flagged."""
    by = defaultdict(lambda: {"signals": [], "chemistry": []})
    for m in moves:
        if m.get("company") and m.get("target"):
            by[(m["company"], m["target"])]["signals"].append(m)
            by[(m["company"], m["target"])]["chemistry"] += m.get("chemistry", [])
    filers = Counter(t for (_, t) in by)                      # distinct companies with a covalent move on each target
    out = []
    for (company, target), v in by.items():
        sigs = sorted(v["signals"], key=lambda s: s["date"] or "9999")
        types = {s["type"] for s in sigs}
        out.append({"company": company, "target": target,
                    "earliest": sigs[0]["date"], "latest": sigs[-1]["date"], "earliest_type": sigs[0]["type"],
                    "cross_linked": len(types) > 1, "has_chemistry": bool(v["chemistry"]),
                    "warheads": sorted(set(v["chemistry"])),
                    "target_filers": filers[target], "lane": _lane(target, filers[target]),
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
    seen = {tuple(x) for x in json.loads(SEEN.read_text())} if SEEN.exists() else set()
    new = {(e["company"], e["target"]) for e in log} - seen
    # Axiom-relevance: lead with WHITE SPACE (open/emerging targets); bury the saturated franchises.
    white = sorted([e for e in log if e["lane"] in ("OPEN", "EMERGING")], key=lambda e: e["latest"], reverse=True)
    contested = [e for e in log if e["lane"] == "CONTESTED"]
    crowded = sorted({e["target"] for e in log if e["lane"] == "CROWDED"})
    clinical = [e for e in log if e["n_clinical"]]

    print(f"=== covalent-move log {today}: {len(log)} moves ===")
    print(f"\n🎯 WHITE SPACE — covalent moves on OPEN/EMERGING targets (the fast-follow gold): {len(white)}")
    for e in white[:22]:
        tag = "🆕 " if (e["company"], e["target"]) in new else "   "
        print(f"  {tag}{e['latest']}  {e['company'][:22]:22} {e['target']:9} [{e['lane']}/{e['target_filers']}] "
              f"{('clin ' + e['signals'][-1]['id']) if e['n_clinical'] else ','.join(e['warheads'])}")
    print(f"\n📊 CONTESTED (4–10 filers, watch): {len(contested)} targets · "
          f"CROWDED/SATURATED (ignore — KRAS/EGFR/BTK class): {len(crowded)} → {crowded[:14]}")
    print(f"\n🧪 CLINICAL covalent signals — tagged by lane (crowded = noise):")
    for e in clinical:
        ncts = [s["id"] for s in e["signals"] if s["type"] == "clinical"][:1]
        print(f"   {e['company'][:22]:22} {e['target']:9} [{e['lane']}]  {ncts}  "
              f"{'<- saturated, low value' if e['lane'] == 'CROWDED' else '<- ACTIONABLE'}")

    SEEN.write_text(json.dumps(sorted([[e["company"], e["target"]] for e in log])))
    json.dump(log, open(HERE / "covalent_move_log.json", "w"), indent=1)
    print(f"\nsaved covalent_move_log.json ({len(log)} moves · {len(white)} white-space)")
    return log


def _gpatent(pn):
    return f"https://patents.google.com/patent/{pn.replace('-', '')}"


def report(n=12, leads=True):
    """Per white-space target, a linked dossier wiring the EXISTING pieces: who's-moving (linked patents/trials),
    a deterministic title-confirms-target check, the lead covalent SMILES (feed_build), and CysDB ligandability."""
    import feed_build
    from cysdb import ligandability
    log = unify()
    feed = {e[0]: e for e in json.load(open(HERE / "covalent_warhead_feed_clean.json"))}
    tg = json.load(open(HERE / "patent_targets.json"))
    cys = ligandability()
    white = sorted([e for e in log if e["lane"] in ("OPEN", "EMERGING")], key=lambda e: e["latest"], reverse=True)[:n]

    out = [f"# Covalent white-space dossier — {datetime.date.today().isoformat()}  ({len(white)} targets)\n"]
    for e in white:
        pats = sorted([s for s in e["signals"] if s["type"] == "patent"], key=lambda x: x["date"], reverse=True)
        lp = pats[0] if pats else None
        title = (feed.get(lp["id"], [None] * 4)[3] if lp else "") or ""
        common = ((tg.get(lp["id"]) or {}).get("common") if lp else None) or e["target"]
        names = e["target"].upper() in title.upper() or common.upper() in title.upper()
        cy = cys.get(e["target"])
        feas = (f"✅ ligandable cysteine ({cy['n_lig_cys']} in CysDB chemoproteomics)" if cy and cy.get("ligandable")
                else "⚠️ profiled in CysDB but NO ligandable cysteine found" if cy
                else "❔ not in CysDB — covalent feasibility unverified")
        lead = feed_build.lead_smiles(lp["id"]) if (leads and lp) else None

        out.append(f"## {e['target']}  —  {e['lane']} · {e['target_filers']} covalent filer(s)")
        out.append(f"- **Who's moving:** {e['company']}"
                   + (f" — [{lp['id']}]({_gpatent(lp['id'])}) ({lp['date']}) · _{title[:64]}_" if lp else ""))
        for s in (s for s in e["signals"] if s["type"] == "clinical"):
            out.append(f"- **Clinical:** [{s['id']}](https://clinicaltrials.gov/study/{s['id']}) ({s['date']})")
        out.append(f"- **Title confirms target:** "
                   + ("✅ yes" if names else f"❌ NO — title doesn't name {e['target']} (likely mis-attribution — verify)"))
        out.append(f"- **Chemistry:** warhead {e['warheads']}"
                   + (f" · lead `{lead}`" if lead
                      else " · lead: US/CN patent (no OPS full-text)" if lp and lp["id"][:2] in ("US", "CN")
                      else " · lead: not resolved"))
        out.append(f"- **Covalent feasibility (CysDB):** {feas}")
        out.append(f"- **Verify:** [gene](https://www.genenames.org/data/gene-symbol-report/#!/symbol/{e['target']}) · "
                   f"[trials](https://clinicaltrials.gov/search?term={e['target']}) · "
                   f"[papers](https://europepmc.org/search?query={e['target']})\n")
    md = "\n".join(out)
    (HERE / f"dossier_{datetime.date.today().isoformat()}.md").write_text(md)
    print(md)
    return md


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
    elif "--report" in sys.argv:
        report()
    else:
        run()
