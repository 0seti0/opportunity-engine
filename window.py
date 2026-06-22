"""Window-assessment for the covalent fast-follow feed. Each 2024+ covalent program is attributed to
its REAL target (SureChEMBL GeneOrProtein annotation, HGNC-resolved -> patent_targets.json), then paired
with two signals to flag whether the best-in-class window is still OPEN or already CLOSED:

  PRIMARY  = covalent-race CROWDING: distinct assignees filing covalent patents on the target NOW
             (current, covalent-specific, no lag -- the PARP7 lesson: a crowded race = closed).
  SECONDARY= CT.gov clinical lead-time: how long the target has been in interventional trials (context;
             lags + biomarker-noisy, so it informs but does not set the headline window).
"""
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

import backtest as bt   # reuse _get, ct_trials, CACHE (cached, retrying)

HERE = Path(__file__).parent
NOW = 2026
LEGAL = {"INC", "LLC", "LTD", "LIMITED", "CO", "CORP", "CORPORATION", "COMPANY", "GMBH", "AG", "SA", "SAS",
         "BV", "NV", "KG", "KGAA", "PLC", "LP", "LLP", "PTE", "PTY", "SRL", "SPA", "OY", "AB", "AS", "ULC"}


# degraders/molecular glues: the warhead match is LINKER / E3-recruiter chemistry, not covalent target
# engagement (verified: IRAK4/STAT6/SMARCA2/MDM2 "acrylamide" hits are all CRBN/VHL degraders). Out of scope
# for a covalent-INHIBITOR radar, so drop them.
DEGRADER = re.compile(r"DEGRADER|PROTAC|MOLECULAR GLUE|BIFUNCTIONAL|CEREBLON|\bCRBN\b|VHL LIGAND|GLUTARIMIDE|\bIMID\b|E3 (UBIQUITIN )?LIGASE")


def norm_assignee(a):                       # one company = one player: drop (US)/(CN) tags + legal suffixes
    a = re.sub(r"\(.*?\)", " ", (a or "").upper())
    a = re.sub(r"[.,/]", " ", a)
    toks = a.split()
    while toks and toks[-1] in LEGAL:
        toks.pop()
    return " ".join(toks).strip()


def earliest_clinical(term):                # earliest year-bracket the target entered interventional trials.
    for hi in (2012, 2015, 2018, 2020, 2022, 2024):
        if bt.ct_trials(term, "2000-01-01", f"{hi}-12-31") > 0:
            return hi
    return None


def run():
    feed = json.load(open(HERE / "covalent_warhead_feed_clean.json"))   # [pn, pub_date, assignees, title, n_warhead, warheads]
    targets = json.load(open(HERE / "patent_targets.json"))             # pn -> {sym, conf, others}

    # a covalent PROGRAM on target X = X named in the patent's title or claims (body-only mentions are
    # too weak: they let polymer/contact-lens false-positives attach to spurious genes).
    recent = defaultdict(list)
    degr = 0
    for e in feed:
        if (e[1] or "")[:4] >= "2024":
            if DEGRADER.search((e[3] or "").upper()):
                degr += 1
                continue
            t = targets.get(e[0])
            if t and t.get("sym") and t.get("conf") in ("title", "claims"):
                recent[t["sym"]].append(e)
    print(f"recent (2024+) covalent-INHIBITOR programs (target in title/claims, degraders excluded): "
          f"{sum(len(v) for v in recent.values())} patents across {len(recent)} targets "
          f"({degr} degrader/glue patents excluded)")

    rows = []
    for i, (sym, es) in enumerate(recent.items()):
        firms = {norm_assignee(e[2].split(';')[0]) for e in es} - {""}   # distinct NORMALIZED firms = the race
        players = len(firms)
        cname = Counter(targets[e[0]].get("common") or sym for e in es).most_common(1)[0][0]  # CT.gov search term
        ey = earliest_clinical(cname)
        rec = bt.ct_trials(cname, "2022-01-01", "2026-12-31")
        clinic_yrs = (NOW - ey) if ey else None
        window = ("CLOSED" if players >= 5 else "CONTESTED" if players >= 3 else "OPEN")   # crowding-primary
        rows.append({"target": sym, "common": cname, "cov_patents": len(es), "players": players,
                     "earliest_trial": ey, "clinic_yrs": clinic_yrs, "recent_trials": rec,
                     "window": window, "latest_pub": max(e[1] for e in es), "assignees": sorted(firms)[:8]})
        if i % 25 == 24:
            bt.CACHE.write_text(json.dumps(bt._cache))
    bt.CACHE.write_text(json.dumps(bt._cache))

    order = {"OPEN": 0, "CONTESTED": 1, "CLOSED": 2}
    rows.sort(key=lambda r: (order[r["window"]], -r["cov_patents"]))
    print("\nPARP7 (the lesson):", next((r for r in rows if r["target"].upper() in ("PARP7", "TIPARP")), "not found"))
    for w in ("OPEN", "CONTESTED", "CLOSED"):
        sub = [r for r in rows if r["window"] == w]
        print(f"\n=== {w}  ({len(sub)} targets) ===")
        for r in sub[:16]:
            cl = f"clinic {r['clinic_yrs']}y" if r["clinic_yrs"] is not None else "no clinical"
            print(f"  {r['target']:12} players={r['players']:2} cov_patents(24+)={r['cov_patents']:3} "
                  f"{cl:12} recent_trials={r['recent_trials']:4}")
    json.dump(rows, open(HERE / "covalent_window.json", "w"), indent=1)
    print(f"\nsaved {len(rows)} target-windows -> covalent_window.json")


if __name__ == "__main__":
    run()
