"""Accurate WHO'S MOVING per target — fixes the patent-feed undercount that made NLRP3 look like
"2 players" when it's a crowded clinical field. Combines:
  clinical_programs() — ClinicalTrials.gov sponsors + drugs + phases (the REAL competition, ALL modalities)
  covalent_filers()   — the covalent-warhead patent feed (covalent-specific IP signal only)
Crowding is judged on the real clinical player count. (Web/pipeline enrichment — e.g. BioAge BGE-102,
which CT.gov keyword search misses — is layered on via claude -p in the report step.)
"""
import json
import re
import urllib.parse
from collections import defaultdict
from pathlib import Path

import backtest as bt

HERE = Path(__file__).parent
PHASE = {"PHASE3": 3, "PHASE2": 2, "PHASE1": 1, "EARLY_PHASE1": 0.5, "NA": 0}
ACADEMIC = ("univers", "hospital", "national", "institut", "college", "medical center", "medical centre",
            "foundation", "health system", "clinic", "cancer center", "assistance publique", "centre hospital")
LEGAL = {"INC", "LLC", "LTD", "LIMITED", "CO", "CORP", "CORPORATION", "COMPANY", "GMBH", "AG", "SA", "SAS",
         "BV", "NV", "KG", "KGAA", "PLC", "LP", "LLP", "PTE", "PTY", "SRL", "SPA", "OY", "AB", "AS", "ULC"}
CN = {"诺华": "NOVARTIS", "吉利德": "GILEAD", "默克": "MERCK", "勃林格殷格翰": "BOEHRINGER INGELHEIM",
      "罗氏": "ROCHE", "辉瑞": "PFIZER", "阿斯利康": "ASTRAZENECA", "基因泰克": "GENENTECH",
      "缬图": "VENTUS", "米拉蒂": "MIRATI", "缆图": "BLUEPRINT", "锐新": "REVOLUTION"}
# collapse a company's many legal-entity spellings to one canonical filer (a keyword -> the canonical name).
# Subsidiary -> ULTIMATE PARENT (point-in-time acquisitions) so one owner isn't counted as several filers —
# that fragmentation inflated crowding lanes and dismissed real white space (audit). Hand-maintained snapshot;
# can only MERGE filers (never split), so it only ever surfaces more white space. The digest also prints the
# distinct filer names at a lane boundary, so a stale mapping is visible to a human.  # ponytail: parent map
CANON = {"GILEAD": "GILEAD SCIENCES", "NOVARTIS": "NOVARTIS", "REVOLUTION": "REVOLUTION MEDICINES",
         "BOEHRINGER": "BOEHRINGER INGELHEIM", "ASTRAZENECA": "ASTRAZENECA", "ROCHE": "ROCHE",
         "GENENTECH": "ROCHE", "PFIZER": "PFIZER", "ARRAY": "PFIZER", "AMGEN": "AMGEN", "JANSSEN": "JANSSEN",
         "BEIGENE": "BEIGENE", "TAKEDA": "TAKEDA", "LILLY": "ELI LILLY", "LOXO": "ELI LILLY",
         "BLUEPRINT": "SANOFI", "INCYTE": "INCYTE", "EISAI": "EISAI", "ABBVIE": "ABBVIE",
         "PHARMACYCLICS": "ABBVIE", "SANOFI": "SANOFI", "VERTEX": "VERTEX", "BAYER": "BAYER", "MERCK": "MERCK",
         "BRISTOL": "BRISTOL MYERS SQUIBB", "MIRATI": "BRISTOL MYERS SQUIBB", "VENTUS": "VENTUS",
         "NIMBUS": "NIMBUS THERAPEUTICS"}                  # per-asset SPVs (Nimbus Lakshmi/Saturn/...) -> one owner
# financial trustees / collateral agents that appear as the assignee when a patent is pledged in debt
# financing — NOT a drug-company filer. Drop them so they neither inflate crowding nor become a phantom mover.
TRUSTEE = re.compile(r"WILMINGTON|ANKURA|GLAS TRUST|COLLATERAL AGENT|ALTER DOMUS|U\.?\s?S\.?\s?BANK|"
                     r"UMB BANK|COMPUTERSHARE|ACQUIOM|WILMINGTON SAVINGS", re.I)


def firm(a):                                       # normalize assignee/sponsor to ONE canonical filer
    a = a or ""
    for cn, en in CN.items():
        a = a.replace(cn, en)                       # CN pharma name -> EN, inline
    parts = [p for p in re.sub(r"[^\x00-\x7F]+", " ", a).split(";") if p.strip()]   # primary assignee, ASCII only
    a = re.sub(r"\(.*?\)|[.,/]", " ", (parts[0] if parts else "").upper())
    if TRUSTEE.search(a):
        return ""                                   # collateral-agent assignee, not a real filer -> callers drop it
    toks = a.split()
    while toks and toks[-1] in LEGAL:               # drop trailing INC/LLC/AG/...
        toks.pop()
    name = " ".join(toks).strip()
    return next((v for k, v in CANON.items() if k in name), name)


def clinical_programs(target):                     # ClinicalTrials.gov industry programs: sponsor, drug, phase
    url = "https://clinicaltrials.gov/api/v2/studies?" + urllib.parse.urlencode({
        "query.term": f"{target} inhibitor", "filter.advanced": "AREA[StudyType]INTERVENTIONAL",
        "pageSize": 150, "format": "json"})
    prog = defaultdict(lambda: {"ph": set(), "drugs": set()})
    for s in bt._get(url).get("studies", []):
        ps = s.get("protocolSection", {})
        spon = ps.get("sponsorCollaboratorsModule", {}).get("leadSponsor", {}).get("name", "")
        if not spon or any(w in spon.lower() for w in ACADEMIC) or " " in spon and spon.istitle() and len(spon.split()) == 2:
            continue                               # drop academic + individual-investigator (Firstname Lastname) sponsors
        prog[spon]["ph"].update(ps.get("designModule", {}).get("phases", []))
        prog[spon]["drugs"].update(i.get("name") for i in ps.get("armsInterventionsModule", {}).get("interventions", [])
                                   if i.get("type") == "DRUG" and "placebo" not in (i.get("name") or "").lower())
    out = []
    for spon, p in prog.items():
        mph = max(p["ph"], key=lambda x: PHASE.get(x, 0)) if p["ph"] else "NA"
        out.append({"sponsor": spon, "phase": mph, "phase_rank": PHASE.get(mph, 0),
                    "drug": (sorted(p["drugs"])[:1] or [None])[0]})
    return sorted(out, key=lambda x: -x["phase_rank"])


def covalent_filers(target):                       # who's filing COVALENT-warhead patents (one slice only)
    from attribute import trusted                   # post-audit grade (title/claims composition, noise removed)
    tg = json.load(open(HERE / "patent_targets.json"))
    feed = {e[0]: e for e in json.load(open(HERE / "covalent_warhead_feed_clean.json"))}
    firms = defaultdict(lambda: {"latest": "", "n": 0})
    for pn, d in tg.items():
        if d.get("sym") == target and trusted(d) and pn in feed:
            e = feed[pn]
            f = firm(e[2].split(';')[0])
            firms[f]["n"] += 1
            firms[f]["latest"] = max(firms[f]["latest"], e[1] or "")
    return sorted(({"firm": f, **v} for f, v in firms.items() if f), key=lambda x: x["latest"], reverse=True)


def web_landscape(target, retries=1):
    """LLM/web leg (claude -p, NO API key) — catches programs CT.gov keyword search misses (e.g. BioAge
    BGE-102) and flags covalent vs non-covalent. Slow (~1 web-search call). [{company,drug,stage,covalent}]."""
    import subprocess

    from llm_claude_code import _extract_json
    prompt = (f"Search the web for {target} inhibitor/modulator drug programs active 2024-2026. Return ONLY a "
              f"JSON array, no prose: [{{\"company\":\"\",\"drug\":\"\",\"stage\":\"\",\"covalent\":false}}]. Be comprehensive.")
    for _ in range(retries + 1):
        try:
            r = subprocess.run(["claude", "-p", prompt, "--output-format", "json"],
                               capture_output=True, text=True, timeout=400)
            out = _extract_json(json.loads(r.stdout).get("result", ""))
            if isinstance(out, list):
                return out
        except Exception:
            pass
    return []


def whos_moving(target, web=False):
    clin = clinical_programs(target)               # ClinicalTrials.gov (deterministic)
    cov = covalent_filers(target)                  # covalent-warhead patents (deterministic)
    land = web_landscape(target) if web else []    # web/pipeline (claude -p) — optional, catches CT.gov gaps
    n = max(len(clin), len(land))
    crowding = "OPEN" if n <= 1 else "CONTESTED" if n <= 4 else "CROWDED"
    cov_clinical = [x for x in land if x.get("covalent")]
    return {"target": target, "crowding": crowding, "n_clinical": len(clin), "n_landscape": len(land),
            "clinical_programs": clin, "covalent_filers": cov, "landscape": land,
            "covalent_lane": "contested" if cov or cov_clinical else "open"}


if __name__ == "__main__":
    import sys
    if "--selftest" in sys.argv:                       # firm() canonicalization incl. subsidiary -> parent
        assert firm("Gilead Sciences, Inc. ; Gilead Sciences Ireland UC") == "GILEAD SCIENCES"
        assert firm("Array Biopharma Inc") == "PFIZER"                       # acquired subsidiary -> parent
        assert firm("Mirati Therapeutics, Inc.") == "BRISTOL MYERS SQUIBB"
        assert firm("Blueprint Medicines Corporation") == "SANOFI"
        assert firm("Genentech, Inc.") == "ROCHE"
        assert firm("吉利德") == "GILEAD SCIENCES"                            # CN -> EN -> canonical
        assert firm("Nimbus Lakshmi, Inc.") == "NIMBUS THERAPEUTICS"         # per-asset SPV -> one owner
        assert firm("Wilmington Trust, National Association") == ""          # financial trustee -> dropped
        assert firm("Ankura Trust Company, LLC") == ""
        assert firm("Some Tiny Biotech LLC") == "SOME TINY BIOTECH"          # unknown -> cleaned passthrough
        print("whos_moving firm() self-check OK")
        sys.exit()
    wm = whos_moving(sys.argv[1] if len(sys.argv) > 1 else "NLRP3")
    print(f"{wm['target']}: {wm['crowding']} — {wm['n_clinical']} industry clinical programs")
    for c in wm["clinical_programs"][:12]:
        print(f"   {c['phase']:8} {c['sponsor'][:34]:34} {c['drug']}")
    print("  covalent IP filers:", [f"{c['firm']}({c['latest'][:7]})" for c in wm["covalent_filers"]])
