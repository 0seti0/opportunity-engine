"""Top up the (lagged) SureChEMBL census with Day-1 hits from the live EPO radar, for a target watchlist.

The bulk census trails reality by SureChEMBL's extraction lag; `feed_build.screen` pulls covalent patents from
EPO OPS AT PUBLICATION. This sweeps a watchlist, runs the radar per target, fetches each hit's biblio, and
merges it into BOTH the feed and patent_targets.json — tagged source='radar', grade='trusted' (feed_build is
the authoritative attributor: the warhead is confirmed on the exemplified composition-of-matter, so these
bypass attribute.py's text-mine entirely). Deduped against the existing census by normalized publication number.

  python census_topup.py                 # sweep watchlist.txt (or the default set)
  python census_topup.py EGFR BTK        # sweep these targets
Cron weekly next to refresh.py: refresh = bulk + lagged; topup = watchlist + Day-1.
ponytail: family-level dedup (WO vs its EP/US member) is by publication number only — a WO hit can still
double a SureChEMBL EP-member of the same family; add family_id dedup if that noise matters.
"""
import json
import re
import sys
import time
from pathlib import Path

import feed_build
import ops

HERE = Path(__file__).parent
CENSUS = HERE / "patent_targets.json"
FEED = HERE / "covalent_warhead_feed_clean.json"
DEFAULT = ["EGFR", "BTK", "KRAS", "BRAF", "FGFR2", "FGFR3", "RET", "JAK3", "SOS1", "CDK7",
           "PRMT5", "WRN", "KEAP1", "NLRP3", "STAT6", "TYK2", "PARP1", "USP1"]


def _norm(pn):
    """pub number -> dedup key (drop dashes + kind): 'EP-4067347-B1' / 'EP4067347A1' -> 'EP4067347'."""
    return re.sub(r"^([A-Z]{2}\d+).*", r"\1", (pn or "").replace("-", ""))


def _merge(census, feed, have, ev, biblio):
    """Add one radar event to (census, feed) if its pub number isn't already present. Returns True if added."""
    if _norm(ev["pn"]) in have:
        return False
    date, appl, title = biblio
    warheads = sorted({w for ld in ev.get("leads", []) for w in ld["warhead"]})
    feed.append([ev["pn"], date, appl, title, ev["n_confirmed"], warheads, "radar"])
    census[ev["pn"]] = {"sym": ev["target"], "common": ev["target"], "conf": "radar",
                        "grade": "trusted", "source": "radar", "warheads": warheads, "others": []}
    have.add(_norm(ev["pn"]))
    return True


def watchlist():
    f = HERE / "watchlist.txt"
    return [g.strip() for g in f.read_text().split() if g.strip()] if f.exists() else DEFAULT


def run(targets=None):
    targets = targets or watchlist()
    census = json.loads(CENSUS.read_text()) if CENSUS.exists() else {}
    feed = json.loads(FEED.read_text()) if FEED.exists() else []
    have = {_norm(pn) for pn in census} | {_norm(e[0]) for e in feed}
    added = 0
    for t in targets:
        for ev in feed_build.screen(t):
            time.sleep(0.5)
            if _merge(census, feed, have, ev, ops.biblio(ev["pn"])):
                added += 1
                e = feed[-1]
                print(f"  + {e[0]:16} {ev['target']:8} {e[1]} {(e[2] or '')[:26]:26} warheads={e[5]}")
    FEED.write_text(json.dumps(feed))
    CENSUS.write_text(json.dumps(census, indent=1))
    print(f"radar top-up: +{added} Day-1 covalent patents (source='radar') merged "
          f"({len(census)} attributed · {len(feed)} feed).")


if __name__ == "__main__":
    if "--selftest" in sys.argv:                # offline check of the merge + dedup
        from attribute import trusted
        c, f, h = {}, [], set()
        ev = {"pn": "WO2026115265", "target": "EGFR", "n_confirmed": 1,
              "leads": [{"name": "x", "smiles": "C=CC(=O)Nc1ccccc1", "warhead": ["acrylamide"]}]}
        assert _merge(c, f, h, ev, ("2026-06-04", "NEOPHORE LTD", "INHIBITOR COMPOUNDS"))
        assert not _merge(c, f, h, {**ev, "pn": "WO-2026115265-A1"}, (None, None, None))   # format variant -> deduped
        assert c["WO2026115265"]["grade"] == "trusted" and c["WO2026115265"]["source"] == "radar"
        assert f[0][5] == ["acrylamide"] and f[0][6] == "radar" and len(f) == 1
        assert trusted(c["WO2026115265"])           # radar entries pass the consumer filter -> flow downstream
        print("census_topup self-check OK:", f[0])
    else:
        run(sys.argv[1:] or None)
