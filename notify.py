"""Weekly Slack digest — conversational covalent-move briefing, web-fact-checked, posted to Slack.

Order:  (1) 👁 YOUR PIPELINE — any new covalent patent / trial on a target in .pipeline (e.g. RAF1), pinned
            FIRST regardless of crowding, because those are your own programs to defend;
        (2) 🆕 new covalent targets I found (sole filer, no prior covalent chemistry);
        (3) other moves;  (4) ⚠️ traps.
Pipeline:  compose -> verify links resolve -> claude -p FACT-CHECK (also WEB-validates each 'only company'
claim) -> re-render softening any overstated claim -> post. Diff-only via .notify_seen.json. No API key
(claude -p on Max). Secrets/state gitignored: .slack_webhook, .pipeline, .notify_seen.json.
"""
import json
import subprocess
import urllib.parse
import urllib.request
from datetime import date, timedelta
from pathlib import Path

HERE = Path(__file__).parent
SEEN = HERE / ".notify_seen.json"
WEBHOOK = HERE / ".slack_webhook"
PIPELINE = HERE / ".pipeline"
MOVES = HERE / "covalent_move_log.json"
INTEL = HERE / ".covalent_intel.json"
WHITE_SPACE = ("OPEN", "EMERGING")


def _gpatent(pn):
    return f"https://patents.google.com/patent/{(pn or '').replace('-', '')}"


def _trial(nct):
    return (f"https://euclinicaltrials.eu/ctis-public/#/study/{nct}" if (nct or '').startswith("20")
            else f"https://clinicaltrials.gov/study/{nct}")


def _link_ok(url):
    try:
        with urllib.request.urlopen(urllib.request.Request(url, method="HEAD"), timeout=20) as r:
            return r.status < 400
    except Exception:
        return False


def _ctgov_recent(target, since):
    """Any interventional trial naming `target`, first-posted since `since` (live CT.gov) — broad, for the
    pipeline watch (we want EVERY movement on our own targets, not just covalent-suspect ones)."""
    q = {"query.term": target, "format": "json", "pageSize": 10,
         "filter.advanced": f"AREA[StudyType]INTERVENTIONAL AND AREA[StudyFirstPostDate]RANGE[{since},MAX]"}
    try:
        with urllib.request.urlopen("https://clinicaltrials.gov/api/v2/studies?" + urllib.parse.urlencode(q), timeout=30) as r:
            d = json.load(r)
    except Exception:
        return []
    out = []
    for s in d.get("studies", []):
        ps = s.get("protocolSection", {}); idm = ps.get("identificationModule", {})
        out.append({"id": idm.get("nctId"), "url": _trial(idm.get("nctId")), "kind": "trial",
                    "company": ps.get("sponsorCollaboratorsModule", {}).get("leadSponsor", {}).get("name", ""),
                    "date": ps.get("statusModule", {}).get("studyFirstPostDateStruct", {}).get("date", "")})
    return out


def pipeline_findings(since, seen):
    """For each target in .pipeline (your internal programs): ANY new TRUSTED covalent patent (census) + ANY
    new trial (live CT.gov), regardless of crowding/lane. These are pinned to the top of the digest."""
    targets = [t.strip().upper() for t in PIPELINE.read_text().split() if t.strip()] if PIPELINE.exists() else []
    if not targets:
        return []
    from attribute import trusted
    tg = json.loads((HERE / "patent_targets.json").read_text())
    feed = {e[0]: e for e in json.loads((HERE / "covalent_warhead_feed_clean.json").read_text())}
    out = []
    for t in targets:
        moves = []
        for pn, d in tg.items():                                    # trusted covalent patents on t (not body-mention noise)
            if (d.get("sym") or "").upper() == t and trusted(d) and (t, pn) not in seen:
                e = feed.get(pn)
                moves.append({"kind": "covalent patent", "id": pn, "url": _gpatent(pn),
                              "company": (e[2] if e else "") or "", "date": (e[1] if e else "") or ""})
        for tr in _ctgov_recent(t, since):                          # any new trial naming t
            if (t, tr["id"]) not in seen:
                moves.append(tr)
        out.append({"target": t, "moves": sorted(moves, key=lambda m: m.get("date", ""), reverse=True)})
    return out


def findings():
    """New (company,target) white-space moves. novel = sole covalent filer AND no covalent chemistry on record."""
    import covindb
    moves = json.loads(MOVES.read_text()) if MOVES.exists() else []
    intel = json.loads(INTEL.read_text()) if INTEL.exists() else {}
    seen = {tuple(x) for x in json.loads(SEEN.read_text())} if SEEN.exists() else set()
    out = []
    for e in moves:
        if (e["company"], e["target"]) in seen or e.get("lane") not in WHITE_SPACE:
            continue
        pats = sorted([s for s in e["signals"] if s["type"] == "patent"], key=lambda s: s["date"], reverse=True)
        clins = sorted([s for s in e["signals"] if s["type"] == "clinical"], key=lambda s: s["date"], reverse=True)
        cy = covindb.lookup(e["target"])
        out.append({"target": e["target"], "company": e["company"],
                    "filers": e.get("target_filer_names", [e["company"]]), "n_filers": e.get("target_filers", 1),
                    "novel": e.get("target_filers", 1) <= 1 and cy is None,
                    "patent": ({"id": pats[0]["id"], "url": _gpatent(pats[0]["id"]), "date": pats[0]["date"]} if pats else None),
                    "trial": ({"id": clins[0]["id"], "url": _trial(clins[0]["id"]), "date": clins[0]["date"]} if clins else None),
                    "date": (pats[0]["date"] if pats else (clins[0]["date"] if clins else "")),
                    "cys": (cy or {}).get("cys_sites", [])[:1], "intel": intel.get(e["target"], {})})
    return out


def _line(f, soften):
    t = f["target"]
    when = (f" — patent {f['patent']['date']}" if f.get("patent") else f" — trial {f['trial']['date']}" if f.get("trial") else "")
    cys = f" It has a reactive cysteine ({f['cys'][0]}), so a covalent approach is feasible." if f["cys"] else ""
    link = (f" <{f['patent']['url']}|see the patent →>" if f.get("patent") and f["patent"]["url"]
            else f" <{f['trial']['url']}|see the trial →>" if f.get("trial") and f["trial"]["url"] else "")
    if t in soften:
        who = f"{f['company']} is filing covalent patents on it (others are too)"
    elif f["n_filers"] <= 1:
        who = f"{f['company']} is the only company filing covalent patents on it in my data"
    else:
        who = f"{f['n_filers']} companies are filing covalent patents on it ({', '.join(f['filers'][:3])})"
    return f"• *{t}*{when} — {who}.{cys}{link}"


def compose(fs, soften=(), pipeline=None):
    """👁 pipeline first, then 🆕 new targets, then other moves, then traps. Every line carries its date."""
    soften = set(soften)
    lines = [f"*🎯 Covalent Radar — {date.today().isoformat()}*",
             "I checked the latest covalent-warhead patents and Phase-1 trials. Here's what's new since my last check (with dates):"]
    if pipeline:
        lines.append("*👁 On your pipeline — watching these closely:*")
        for p in pipeline:
            if p["moves"]:
                for m in p["moves"][:4]:
                    lines.append(f"• *{p['target']}* — new {m['kind']} from {(m.get('company') or '?')[:26]} "
                                 f"({m.get('date', '')}). <{m['url']}|see it →>")
            else:
                lines.append(f"• *{p['target']}* — no new covalent patents or trials since last check (still clear).")
    traps = [f for f in fs if f["intel"].get("verdict") == "TRAP"]
    novel = [f for f in fs if f["novel"] and f["target"] not in soften and f not in traps]
    others = [f for f in fs if f not in novel and f not in traps]
    if novel:
        lines.append("*🆕 New covalent targets I found* — sole covalent filer, and no prior covalent chemistry on record:")
        lines += [_line(f, soften) for f in novel]
    if others:
        lines.append("*Other covalent moves:*")
        lines += [_line(f, soften) for f in others]
    for f in traps:
        adv = f["intel"].get("most_advanced", "an existing program")
        res = " and on-target resistance has been reported" if f["intel"].get("resistance_reported") else ""
        when = (f"patent {f['patent']['date']}" if f.get("patent") else f"trial {f['trial']['date']}" if f.get("trial") else "")
        link = (f" <{f['patent']['url']}|see the patent →>" if f.get("patent") and f["patent"]["url"]
                else f" <{f['trial']['url']}|see the trial →>" if f.get("trial") and f["trial"]["url"] else "")
        lines.append(f"⚠️ *{f['target']}* ({when}) — {f['company']} has a new covalent move, but I'd skip it: a "
                     f"covalent {f['target']} inhibitor is already clinical ({adv}){res}, so the angle is taken.{link}")
    return "\n\n".join(lines)


def web_validate(sole_targets):
    """The LLM fact-check before sending. The composed message is deterministic from our data (it can't
    misstate a fact), so the one claim that needs checking against REALITY is 'X is the only company filing
    covalent patents on it' — our census is Western-skewed and lagged, so it can be wrong. claude -p searches
    the web for each and returns the targets to SOFTEN (where other covalent filers exist). Small, focused
    prompt -> fast and reliable. On ANY failure, softens ALL of them (never post an unverified 'only' claim)."""
    if not sole_targets:
        return {"soften": [], "checked": True}
    from llm_claude_code import _extract_json
    prompt = ("You know the covalent drug-discovery landscape well. For each target below, do you KNOW that two "
              "or more different companies already have a covalent (irreversible-warhead) small-molecule program "
              f"(patent or clinical) against it? Targets: {sole_targets}. Only include a target if you are "
              "CONFIDENT it has multiple covalent players; leave out ones you're unsure about (we keep an 'in my "
              "data' qualifier for those). You may do a quick web check, but answer from knowledge if you can.\n"
              'Return ONLY JSON: {"multiple_filers": ["TARGET", ...]}.')
    import re
    try:
        r = subprocess.run(["claude", "-p", prompt, "--output-format", "json"], capture_output=True, text=True, timeout=300)
        res = json.loads(r.stdout).get("result", "")          # the model often appends prose after the JSON,
        m = re.search(r'multiple_filers"?\s*:\s*\[([^\]]*)\]', res)   # so pull the array directly (prose-robust)
        if m is not None:
            flagged = re.findall(r'"([A-Za-z0-9]+)"', m.group(1))
            return {"soften": [t for t in flagged if t in sole_targets], "checked": True}
    except Exception:
        pass
    return {"soften": list(sole_targets), "checked": False}


def post_slack(text):
    req = urllib.request.Request(WEBHOOK.read_text().strip(), data=json.dumps({"text": text}).encode(),
                                 headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=30) as r:
        return r.status


def run(dry=False, limit=10):
    since = (date.today() - timedelta(days=90)).isoformat()
    seen = {tuple(x) for x in json.loads(SEEN.read_text())} if SEEN.exists() else set()
    pipeline = pipeline_findings(since, seen)
    fs = findings()
    if not fs and not any(p["moves"] for p in pipeline):
        print("Nothing new since last check."
              + (" Pipeline still clear: " + ", ".join(p["target"] for p in pipeline) if pipeline else ""))
        return
    fs.sort(key=lambda f: (f["intel"].get("verdict") == "TRAP", not f["novel"], f["date"] < "0"))
    fs.sort(key=lambda f: f["date"], reverse=True)
    fs.sort(key=lambda f: (f["intel"].get("verdict") == "TRAP", not f["novel"]))   # novel block, others, traps
    shown, overflow = fs[:limit], max(0, len(fs) - limit)
    for f in shown:
        for kind in ("patent", "trial"):
            if f.get(kind) and not _link_ok(f[kind]["url"]):
                f[kind]["url"] = ""
    val = web_validate([f["target"] for f in shown if f["n_filers"] <= 1])   # LLM fact-check of the 'only company' claims
    final = compose(shown, soften=val["soften"], pipeline=pipeline)
    if overflow:
        final += f"\n\n_(+{overflow} more new moves I'm tracking — I'll surface them as they develop.)_"
    if not val["checked"]:
        final += "\n\n_(Couldn't web-verify the sole-filer claims this run, so I've softened them to be safe.)_"
        print("(web check unavailable -> softened all sole-filer claims)")
    elif val["soften"]:
        print("(softened after web check:", val["soften"], ")")
    if dry or not WEBHOOK.exists():
        print(("=== DRY RUN" + (" — no .slack_webhook, not posting" if not WEBHOOK.exists() else "") + " ===\n") + final)
    else:
        print("posted to Slack (status", post_slack(final), ")")
    seen |= {(f["company"], f["target"]) for f in shown}
    seen |= {(p["target"], m["id"]) for p in pipeline for m in p["moves"]}          # don't re-alert a pipeline move
    SEEN.write_text(json.dumps(sorted(list(k) for k in seen)))


if __name__ == "__main__":
    import sys
    if "--selftest" in sys.argv:
        fs = [{"target": "SARM1", "company": "ELI LILLY", "filers": ["ELI LILLY"], "n_filers": 1, "novel": True,
               "patent": {"id": "WO-2023177972-A1", "url": "https://patents.google.com/patent/WO2023177972A1", "date": "2023-09-21"},
               "trial": None, "date": "2023-09-21", "cys": [], "intel": {}},
              {"target": "WRN", "company": "GILEAD SCIENCES", "filers": ["GILEAD SCIENCES", "INCYTE"], "n_filers": 2,
               "novel": False, "patent": {"id": "US-20250230168-A1", "url": "https://patents.google.com/patent/US20250230168A1", "date": "2025-07-17"},
               "trial": None, "date": "2025-07-17", "cys": [], "intel": {"verdict": "TRAP", "most_advanced": "VVD-214 (Vividion)", "resistance_reported": True}}]
        pipe = [{"target": "RAF1", "moves": []},
                {"target": "BRAF", "moves": [{"kind": "trial", "id": "NCT09", "url": "https://clinicaltrials.gov/study/NCT09", "company": "PFIZER", "date": "2026-06-01"}]}]
        m = compose(fs, pipeline=pipe)
        head = m.split("*🆕")[0]
        assert "On your pipeline" in head and "RAF1" in head and "still clear" in head    # RAF1 pinned first, even with no moves
        assert "BRAF" in head and "NCT09" in head                                          # a pipeline move is shown
        assert m.index("RAF1") < m.index("SARM1")                                          # pipeline before new targets
        assert "ELI LILLY is the only company" in m and "2023-09-21" in m and "VVD-214" in m
        assert "OPEN" not in m and "EMERGING" not in m
        sline = next(l for l in compose(fs, soften={"SARM1"}, pipeline=pipe).split("\n\n") if l.startswith("• *SARM1*"))
        assert "only company" not in sline and "others are too" in sline
        print("notify self-check OK\n" + "-" * 60 + "\n" + m)
        sys.exit()
    run(dry="--dry" in sys.argv)
