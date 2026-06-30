"""Leg 2 — clinical early-warning radar. Catches covalent moves BEFORE the patent publishes, from clinical
registries (no chemistry). New Phase-1 industry trials -> classify covalency + score TARGET NOVELTY.
  (Step 1)  diff new trials — ClinicalTrials.gov + EU CTIS  (WHO ICTRP / ChiCTR have no public API)
  (Step 2a) patent cross-link to our covalent census   (2b) web check (claude -p)   (2c) fingerprint / pubs
  (novelty) cross-ref census + CovalentInDB + prior CT.gov trials -> NOVEL / EMERGING / KNOWN target
  (conf)    earliest signal: AACR/ASCO/ESMO + company-PR sweep for NOVEL targets (claude -p, citation-guarded)

  python clinical_radar.py            # CT.gov + CTIS, cross-link + novelty (fast)
  python clinical_radar.py --web      # + code-name web check (claude -p per drug)
  python clinical_radar.py --conf     # + conference/PR sweep on NOVEL targets (claude -p)
  python clinical_radar.py --selftest # offline checks
"""
import json
import re
import sys
import urllib.parse
import urllib.request
from collections import defaultdict
from datetime import date, timedelta
from pathlib import Path

import backtest as bt
from attribute import trusted
from identity import build_index, load_rows, resolve
from tractable import covalently_tractable
from whos_moving import firm

HERE = Path(__file__).parent
SEEN = HERE / ".clinical_seen.json"          # NCTs already reported (the diff state)
CTGOV = "https://clinicaltrials.gov/api/v2/studies"
HGNC_URL = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
# uppercase tokens that look like genes but are clinical/disease/method acronyms (the attribution lesson)
CLIN_STOP = {"PHASE", "STUDY", "DRUG", "ORAL", "BID", "QD", "FDA", "MTD", "DLT", "PK", "PD", "AE", "SAE",
             "DOSE", "ARM", "USA", "MRI", "PET", "ECOG", "NSCLC", "SCLC", "AML", "CLL", "CML", "HCC", "TNBC",
             "DNA", "RNA", "CSF", "ITT", "OPEN", "LABEL", "SOLID", "ADULT", "ADULTS", "SAFETY", "ADME", "IND",
             "DMD", "ALS", "COPD", "NASH", "MASH", "GVHD", "HBV", "HIV"}
# a covalent SMALL MOLECULE can't be any of these modalities -> drop the trial (intervention name + title)
_NONSM = re.compile(r"\b\d{2,3}Lu\b|\b\d{2,3}Ga\b|\b\d{3}Ac\b|177Lu|68Ga|225Ac|\bLu-|\bGa-|PSMA|radioligand|"
                    r"radiopharma|\bCAR[- ]?T\b|chimeric antigen|lymphodeplet|fludarabine|mab\b|bispecific|"
                    r"\bantibody\b|conjugate|\bADC\b|siRNA|antisense|oligonucleotide|\bmRNA\b|vaccine|"
                    r"cell therapy|\bTCR\b|CRISPR|peptide", re.I)
_PROBE = re.compile(r"bioavailab|drug.drug interaction|\bDDI\b|mass balance|food effect|midazolam|cocktail",
                    re.I)
# targets/mutations that are covalent BY DEFINITION (the warhead is class-required) -> flag covalency from the
# TARGET, bypassing the brand-new-code-name problem (HS-10541, KQB368 — web can't confirm an undisclosed warhead)
_COVALENT_CLASS = re.compile(r"\bG12C\b|\bG12S\b|exon\s?20\s?ins", re.I)


def _smallmol_program(tr, targets):
    """Keep covalently-tractable targets of a SMALL-MOLECULE trial. Modality is judged PER DRUG (a small
    molecule + antibody combo is kept on its small-molecule arm), not nuked because 'antibody' is in the title."""
    if _PROBE.search(" ".join(tr["drugs"] + [tr["title"]])):
        return []
    if not tr["drugs"] or all(_NONSM.search(d or "") for d in tr["drugs"]):   # every drug is a non-SM modality
        return []
    return [t for t in targets if covalently_tractable(t)]


def new_trials(since=None, max_pages=12):
    """ALL Phase-1/Early-Phase-1 industry interventional trials first-posted since `since` (default 90 days),
    PAGINATED — not just the most-recent 200 (that window missed the late-May covalent G12C/BTK starts)."""
    since = since or (date.today() - timedelta(days=90)).isoformat()
    base = ("AREA[Phase](PHASE1 OR EARLY_PHASE1) AND AREA[StudyType]INTERVENTIONAL AND "
            f"AREA[LeadSponsorClass]INDUSTRY AND AREA[StudyFirstPostDate]RANGE[{since},MAX]")
    out, token = [], None
    for _ in range(max_pages):
        q = {"filter.advanced": base, "pageSize": 100, "format": "json"}
        if token:
            q["pageToken"] = token
        r = bt._get(f"{CTGOV}?{urllib.parse.urlencode(q)}")
        for s in r.get("studies", []):
            ps = s.get("protocolSection", {}); idm = ps.get("identificationModule", {})
            out.append({"nct": idm.get("nctId"),
                        "posted": ps.get("statusModule", {}).get("studyFirstPostDateStruct", {}).get("date", ""),
                        "sponsor": ps.get("sponsorCollaboratorsModule", {}).get("leadSponsor", {}).get("name", ""),
                        "title": idm.get("briefTitle", ""),
                        "conditions": " ".join(ps.get("conditionsModule", {}).get("conditions", [])),
                        "drugs": [i.get("name") for i in ps.get("armsInterventionsModule", {}).get("interventions", [])
                                  if i.get("type") == "DRUG" and "placebo" not in (i.get("name") or "").lower()],
                        "summary": (ps.get("descriptionModule", {}).get("briefSummary", "") or "")[:800],
                        "design_text": _design_text(ps)})
        token = r.get("nextPageToken")
        if not token:
            break
    return out


# --- (a) ex-US coverage. EU CTIS has a real public JSON API (verified). WHO ICTRP = no public API
# (registration-gated weekly XML export); ChiCTR = anti-bot HTML — both left as documented slots, NOT faked.
# CTIS gives the EU Phase-1 covalent starts CT.gov misses, mapped to the SAME trial dict so the cascade is reused.
CTIS = "https://euclinicaltrials.eu/ctis-public-api/search"


def _ctis_post(body):
    try:
        req = urllib.request.Request(CTIS, data=json.dumps(body).encode(),
                                     headers={"Content-Type": "application/json"})
        with urllib.request.urlopen(req, timeout=60) as r:
            return json.load(r)
    except Exception:
        return {}                                            # CTIS outage must not break the CT.gov leg


def _ctis_date(s):
    """CTIS dates are 'DD/MM/YYYY' (sometimes prefixed 'GR: ') -> 'YYYY-MM-DD', or '' if unparseable."""
    m = re.search(r"(\d{2})/(\d{2})/(\d{4})", s or "")
    return f"{m.group(3)}-{m.group(2)}-{m.group(1)}" if m else ""


def new_trials_ctis(since=None, max_pages=6):
    """EU CTIS Phase-1 inhibitor trials decided since `since` — EU coverage CT.gov misses, in new_trials()'s
    dict shape so targets_in + the covalency cascade run unchanged. CTIS search needs a text criterion;
    'inhibitor' is the broad net (mirrors the OPS ti=inhibitor precision/recall tradeoff)."""
    since = since or (date.today() - timedelta(days=90)).isoformat()
    out = []
    for page in range(1, max_pages + 1):
        data = _ctis_post({"pagination": {"page": page, "size": 50},
                           "sort": {"property": "decisionDate", "direction": "DESC"},
                           "searchCriteria": {"containAll": "inhibitor"}}).get("data") or []
        if not data:
            break
        for t in data:
            posted = _ctis_date(t.get("decisionDateOverall") or t.get("decisionDate"))
            if posted and posted < since:
                continue                                     # out of window (no early-break: sort field != filter field)
            if not re.search(r"\bphase i\b", t.get("trialPhase", "") or "", re.I):
                continue                                     # Phase-1 only ('Phase I/II' yes; 'Phase II/III/IV' no — \b, not substring)
            if "pharmaceutical" not in (t.get("sponsorType", "") or "").lower():
                continue                                     # industry only (matches CT.gov LeadSponsorClass=INDUSTRY)
            drugs = [p.strip() for p in (t.get("product", "") or "").split(",")
                     if p.strip() and "placebo" not in p.lower()]
            out.append({"nct": t.get("ctNumber", ""), "posted": posted, "sponsor": t.get("sponsor", ""),
                        "title": t.get("ctTitle", ""),
                        "conditions": " ".join(filter(None, [t.get("conditions", ""), str(t.get("therapeuticAreas", ""))])),
                        "drugs": drugs, "summary": (t.get("primaryEndPoint", "") or "")[:800],
                        "design_text": " ".join(filter(None, [t.get("primaryEndPoint", ""), t.get("endPoint", "")]))})
    return out


# --- (b) novelty: how NEW is a covalent move on a target? cross-ref the three prior-art signals we already hold.
def _prior_trials(target, since):
    """CT.gov count of INTERVENTIONAL trials on `target` first-posted BEFORE `since` — its clinical history."""
    q = {"query.term": target, "countTotal": "true", "pageSize": 1, "format": "json",
         "filter.advanced": f"AREA[StudyType]INTERVENTIONAL AND AREA[StudyFirstPostDate]RANGE[MIN,{since}]"}
    return bt._get(f"{CTGOV}?{urllib.parse.urlencode(q)}").get("totalCount", 0)


def _novelty_tier(ip, chem, prior):
    return ["NOVEL", "EMERGING", "KNOWN", "KNOWN"][sum([bool(ip), bool(chem), prior > 0])]


def novelty(target, covidx, since):
    """{covalent_ip (census), covalent_chem (CovalentInDB), prior_trials (CT.gov), tier}. NOVEL = no covalent
    IP + no known covalent chemistry + no prior clinical history = a genuinely new covalent target (the prize)."""
    from covindb import lookup as covindb_lookup
    ip = bool(covidx.get(target))
    chem = covindb_lookup(target) is not None
    prior = _prior_trials(target, since)
    return {"covalent_ip": ip, "covalent_chem": chem, "prior_trials": prior, "tier": _novelty_tier(ip, chem, prior)}


def _design_text(ps):
    """Outcome measures (kept WHOLE — short, high-signal, often 15+ of them) + a capped detailed description.
    Truncating the combined blob dropped late outcomes like a secondary BTK-occupancy measure -> false negatives."""
    oc = ps.get("outcomesModule", {})
    outcomes = " ".join(f"{o.get('measure','')} {o.get('description','') or ''} {o.get('timeFrame','') or ''}"
                        for o in oc.get("primaryOutcomes", []) + oc.get("secondaryOutcomes", []))
    return outcomes + " " + (ps.get("descriptionModule", {}).get("detailedDescription", "") or "")[:4000]


def trial_design(nct):
    """Design text for one trial by NCT (single-study CT.gov endpoint) — for the fingerprint / audit."""
    return _design_text(bt._get(f"{CTGOV}/{nct}?format=json").get("protocolSection", {}))


# 2c — clinical-DESIGN fingerprint: covalency tells that need NO chemistry. The COMBINATION is specific;
# any one alone is weak (occupancy assays aren't covalent-exclusive).
# AUDIT (24 labeled Ph1 trials, adversarially verified): on the CT.gov REGISTRY this fires on only 1/12
# covalent drugs (8% recall, 100% precision) — the tells are almost never in the registry; they live in the
# protocol/publications (10/12 covalent vs 5/12 non-covalent have them there). So 2c-registry is a rare
# BONUS, not a primary detector: lean on 2a (patent cross-link) + 2b (web). A high-recall 2c would fingerprint
# FULL-TEXT publications (PMID -> Europe PMC) and require the covalent-SPECIFIC tells (sustained-PD-past-PK /
# mass-spec adduct / explicit irreversible), since plain occupancy also shows up for reversible drugs.
_FINGERPRINT = {
    "occupancy_timecourse": (3, r"target occupancy|receptor occupancy|target engagement|\boccupancy\b"),
    "sustained_pd":         (2, r"resynthesis|duration of (?:target )?(?:inhibition|engagement)|"
                                r"recovery of \w+ (?:activity|protein|function)|"
                                r"sustained (?:target )?(?:inhibition|pharmacodynamic|engagement)|return to baseline"),
    "adduct_massspec":      (3, r"covalent (?:bond|adduct|occupancy)|activity-based protein profiling|\bABPP\b|"
                                r"probe[- ]?(?:based )?occupancy|mass spectrometr\w+[^.]{0,40}(?:occupancy|adduct|engagement)"),
    "pbmc_pd":              (1, r"\bPBMCs?\b"),
    "covalent_stated":      (3, r"irreversibl\w+|covalent\w*"),
}
_FP = {k: (w, re.compile(p, re.I)) for k, (w, p) in _FINGERPRINT.items()}
FP_THRESHOLD = 4


def fingerprint(design_text):
    """-> {score, signals}. Covalent-suspect when score >= FP_THRESHOLD."""
    hits = {k: w for k, (w, rx) in _FP.items() if rx.search(design_text or "")}
    return {"score": sum(hits.values()), "signals": sorted(hits)}


# full-text-2c — covalency from the trial's PUBLICATIONS (Europe PMC abstracts), since the audit showed the
# tells live in pubs not the registry: this lifts 2c recall 8% -> ~58% at ~88% precision (24-trial audit).
# Negation-guarded so background/comparator mentions ("non-covalent", "previous covalent BTKi") don't fire.
EUPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest"
_PUB_NEG = re.compile(r"non-?covalent|\breversible\b", re.I)
_PUB_STRONG = re.compile(r"covalent(?:ly)? (?:bind|bond|bound|adduct|modif)|\bCys(?:teine)?[\s-]?\d{2,4}\b|"
                         r"\bwarhead\b|acrylamide|Michael acceptor|butynamide|fluoroacrylamide|cyanoacryl", re.I)
_PUB_WORD = re.compile(r"irreversibl\w+|(?<!non[- ])(?<!non)covalent(?:ly)?\b", re.I)


def _pub_verdict(text):
    """A hard covalent MECHANISM phrase fires regardless; a bare 'irreversible/covalent' word fires only if
    the text doesn't also describe the drug as non-covalent/reversible (kills comparator false positives)."""
    if _PUB_STRONG.search(text):
        return True
    return bool(_PUB_WORD.search(text) and not _PUB_NEG.search(text))


def pub_covalent(nct):
    """full-text-2c: is the trial's drug covalent per its linked PUBLICATIONS? Deterministic, no LLM."""
    refs = bt._get(f"{CTGOV}/{nct}?format=json").get("protocolSection", {}).get("referencesModule", {}).get("references", [])
    text = ""
    for pmid in [r["pmid"] for r in refs if r.get("pmid")][:5]:
        res = bt._get(f"{EUPMC}/search?query=EXT_ID:{pmid}%20AND%20SRC:MED&resultType=core&format=json"
                      ).get("resultList", {}).get("result", [])
        if res:
            text += " " + res[0].get("title", "") + " " + (res[0].get("abstractText", "") or "")
    m = _PUB_STRONG.search(text) or _PUB_WORD.search(text)
    return {"covalent": _pub_verdict(text), "evidence": (m.group(0)[:40] if m else "")}


def gene_index():
    """HGNC index (download-if-missing) -> (Index, id->official-symbol)."""
    p = Path("/tmp/hgnc.tsv")
    if not p.exists():
        print("  downloading HGNC complete set ...", flush=True)
        urllib.request.urlretrieve(HGNC_URL, str(p))
    idx = build_index(load_rows(str(p)))
    return idx, {v: k for k, v in idx.symbols.items()}


def targets_in(trial, idx, id2sym):
    """Canonical gene symbols named in the trial (uppercase tokens resolving to a unique HGNC gene)."""
    text = " ".join([trial["title"], trial.get("conditions", "")] + trial["drugs"] + [trial["summary"]])
    hits = set()
    for m in re.finditer(r"\b[A-Z][A-Z0-9]{2,7}\b", text):
        tok = m.group(0)
        if tok in CLIN_STOP:
            continue
        if re.match(r"\s*[- ]\s*(?:neg|pos|low|high|mut|wild|wt\b|amplif|alter|defic)",
                    text[m.end():m.end() + 12], re.I):
            continue                         # biomarker STATUS ("HER2-negative") — a selection criterion, not the target
        r = resolve(tok, idx)
        if r.hgnc_id:                        # RESOLVED (unique); AMBIGUOUS collisions auto-skipped
            hits.add(id2sym[r.hgnc_id])
    return sorted(hits)


def covalent_index():
    """target symbol -> set of firms holding a TRUSTED covalent patent on it (incl. the radar top-up)."""
    tg = json.load(open(HERE / "patent_targets.json"))
    feed = {e[0]: e for e in json.load(open(HERE / "covalent_warhead_feed_clean.json"))}
    idx = defaultdict(set)
    for pn, d in tg.items():
        if d.get("sym") and trusted(d):
            e = feed.get(pn)
            idx[d["sym"]].add(firm(e[2]) if e and e[2] else "")
    return idx


def cross_link(trial, targets, covidx):
    """2a: does the trial's sponsor (or anyone) hold covalent IP on a target it names?"""
    sp = firm(trial["sponsor"])
    for t in targets:
        if sp and sp in covidx.get(t, set()):
            return {"level": "LIKELY", "why": f"sponsor {sp} holds covalent IP on {t}", "target": t}
    for t in targets:
        if covidx.get(t):
            return {"level": "POSSIBLE", "why": f"covalent IP on {t} exists (other filers)", "target": t}
    return {"level": "—", "why": "no covalent patent cross-link", "target": targets[0] if targets else None}


def web_covalent(drug):
    """2b: live web check on the code-name (claude -p, no API key). Code-names aren't in LLM training data,
    so this MUST search the web rather than recall."""
    import subprocess

    from llm_claude_code import _extract_json
    prompt = (f"Is the investigational drug '{drug}' a COVALENT / irreversible inhibitor or a covalent "
              f"degrader? Search the web. Return ONLY JSON: {{\"covalent\": true, \"target\": \"\", "
              f"\"evidence\": \"\", \"sources\": [\"url\"]}}  (covalent null if you can't tell; sources = the "
              f"URLs you actually opened, [] if you opened none).")
    try:
        r = subprocess.run(["claude", "-p", prompt, "--output-format", "json"],
                           capture_output=True, text=True, timeout=300)
        o = _extract_json(json.loads(r.stdout).get("result", ""))
        return o if isinstance(o, dict) else {}
    except Exception:
        return {}


# --- (c) conference/PR sweep: the EARLIEST signal — covalent programs announced at AACR/ASCO/ESMO or in
# company PR before they hit any trial registry. claude -p web, citation-guarded, cached. Best-effort, opt-in.
_CONF_CACHE = HERE / ".conference_pr.json"


def conference_pr(target):
    """Recent (~12mo) conference/PR covalent program on `target` not yet in the registries. claude -p web,
    sources REQUIRED (a 'found' with no source is downgraded), cached. No API key."""
    cache = json.loads(_CONF_CACHE.read_text()) if _CONF_CACHE.exists() else {}
    if target in cache:
        return cache[target]
    import subprocess

    from llm_claude_code import _extract_json
    prompt = (f"Search the web (last 12 months). Is there a COVALENT / irreversible small-molecule program "
              f"against {target} announced in a CONFERENCE abstract (AACR/ASCO/ESMO/AACR-NCI-EORTC) or a company "
              f"press-release / pipeline page, but NOT yet registered on ClinicalTrials.gov or EU CTIS? Return "
              f"ONLY JSON: {{\"found\": true, \"stage\": \"preclinical|IND-enabling|conference\", \"company\": \"\", "
              f"\"event\": \"\", \"evidence\": \"\", \"sources\": [\"url\"]}}. sources = URLs you actually opened "
              f"([] if none); found=false if nothing covalent-specific.")
    try:
        r = subprocess.run(["claude", "-p", prompt, "--output-format", "json"],
                           capture_output=True, text=True, timeout=300)
        o = _extract_json(json.loads(r.stdout).get("result", ""))
        if isinstance(o, dict):
            if o.get("found") and not o.get("sources"):              # citation guard: uncited != found
                o = {**o, "found": False, "evidence": "unverified — no source cited"}
            cache[target] = o
            _CONF_CACHE.write_text(json.dumps(cache, indent=1))
            return o
    except Exception:
        pass
    return {"found": False, "evidence": "sweep unavailable"}


def _dedup_moves(events):
    """Cross-registry dedup: the same covalent move can arrive as an NCT (CT.gov) AND a CTIS number — collapse
    same-(firm, target) events to one, keeping the strongest covalency then the EARLIEST date (the earliest
    signal is the point). Also folds same-(firm,target) intra-registry duplicates (mirrors unify's aggregation)."""
    rank = {"LIKELY": 0, "FINGERPRINT": 1, "POSSIBLE": 2, "—": 3}
    best = {}
    for e in events:
        k = (firm(e.get("sponsor", "")), e.get("target"))
        score = (rank.get(e["covalency"], 3), e.get("posted") or "9999")
        if k not in best or score < (rank.get(best[k]["covalency"], 3), best[k].get("posted") or "9999"):
            best[k] = e
    return list(best.values())


def run(web=False, pubs=False, conf=False):
    since = (date.today() - timedelta(days=90)).isoformat()
    seen = set(json.loads(SEEN.read_text())) if SEEN.exists() else set()
    idx, id2sym = gene_index()
    covidx = covalent_index()
    trials = new_trials(since) + new_trials_ctis(since)      # (a) CT.gov + EU CTIS (ICTRP/ChiCTR: no public API)
    events = []
    for tr in trials:
        if tr["nct"] in seen or not tr["drugs"]:
            continue                         # seen, or a diagnostic/imaging/device trial (no DRUG intervention)
        targets = _smallmol_program(tr, targets_in(tr, idx, id2sym))
        if not targets:
            continue                         # no covalently-tractable small-molecule target -> 2a/2b can't classify
        xl = cross_link(tr, targets, covidx)
        ev = {**tr, "targets": targets, "covalency": xl["level"], "why": xl["why"], "target": xl["target"]}
        if _COVALENT_CLASS.search(" ".join([tr["title"], tr.get("conditions", "")] + tr["drugs"])):
            ev["covalency"] = "LIKELY"                       # covalent-by-class (e.g. all KRAS G12C inhibitors are covalent)
            ev["why"] = "covalent-by-class (G12C/G12S/exon20ins — covalent-definitional mechanism)"
        fp = fingerprint(tr.get("design_text", ""))          # 2c-registry: design tells (cheap; ~8% recall)
        ev["fingerprint"] = fp
        if fp["score"] >= FP_THRESHOLD and ev["covalency"] in ("—", "POSSIBLE"):
            ev["covalency"] = "FINGERPRINT"
            ev["why"] = f"design fingerprint (score {fp['score']}: {', '.join(fp['signals'])})"
        if pubs and ev["covalency"] in ("—", "POSSIBLE", "FINGERPRINT"):   # full-text-2c (deterministic, ~58% recall)
            pc = pub_covalent(tr["nct"])
            if pc["covalent"]:
                ev["covalency"], ev["why"] = "LIKELY", f"pubs: covalent ('{pc['evidence']}')"
        if web and ev["covalency"] != "LIKELY" and tr["drugs"]:           # 2b web (LLM; most accurate, context-aware)
            w = web_covalent(tr["drugs"][0])
            if w.get("covalent") is True and w.get("sources"):            # only a WEB-CITED verdict earns LIKELY
                ev["covalency"], ev["why"] = "LIKELY", f"web: {(w.get('evidence') or '')[:70]} [{len(w['sources'])} src]"
            elif w.get("covalent") is True and ev["covalency"] == "—":     # covalent claimed but uncited -> can't trust
                ev["covalency"], ev["why"] = "POSSIBLE", "web: covalent claimed, NO source cited (unverified)"
        ev.pop("design_text", None)                          # drop the bulky text from the event
        if ev["covalency"] != "—":
            ev["novelty"] = novelty(ev["target"], covidx, since)         # (b) how NEW is this covalent target?
            if conf and ev["novelty"]["tier"] == "NOVEL":               # (c) earliest-signal sweep, only for the prizes
                ev["conference"] = conference_pr(ev["target"])
        events.append(ev)

    n_raw = len(events)
    events = _dedup_moves(events)                            # collapse CT.gov+CTIS duplicates of one program
    nov_rank = {"NOVEL": 0, "EMERGING": 1, "KNOWN": 2}
    rank = {"LIKELY": 0, "FINGERPRINT": 1, "POSSIBLE": 2, "—": 3}
    events.sort(key=lambda e: (rank.get(e["covalency"], 3), nov_rank.get((e.get("novelty") or {}).get("tier"), 3)))
    for e in events:
        if e["covalency"] != "—":
            nov = (e.get("novelty") or {}).get("tier", "?")
            print(f"  [{e['covalency']:8}|{nov:8}] {(e['nct'] or '?')[:18]:18} {e['posted']} {firm(e['sponsor'])[:18]:18} "
                  f"{(e['drugs'][0] if e['drugs'] else '?')[:20]:20} {e['target']}  ({e['why']})")
    SEEN.write_text(json.dumps(sorted(seen | {t['nct'] for t in trials if t['nct']})))
    json.dump(events, open(HERE / "clinical_events.json", "w"), indent=1)
    n_cov = sum(1 for e in events if e["covalency"] != "—")
    n_novel = sum(1 for e in events if e["covalency"] != "—" and (e.get("novelty") or {}).get("tier") == "NOVEL")
    n_ctis = sum(1 for t in trials if (t.get("nct") or "").startswith("20"))     # CTIS ids look like '2025-...'
    print(f"\nclinical radar: {len(trials)} new Ph1 trials ({n_ctis} via EU CTIS) · {len(events)} unique "
          f"(company,target) moves ({n_raw} before cross-registry dedup) · {n_cov} covalent-suspect "
          f"({n_novel} NOVEL targets).")
    return events


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        from identity import Index
        idx = Index(symbols={"EGFR": "HGNC:3236", "BTK": "HGNC:1133", "KRAS": "HGNC:6407"},
                    aliases={"HER2": {"HGNC:3430"}})
        id2 = {"HGNC:3236": "EGFR", "HGNC:1133": "BTK", "HGNC:6407": "KRAS", "HGNC:3430": "ERBB2"}
        tr = {"title": "A Phase 1 Study of an EGFR Inhibitor in NSCLC", "drugs": ["ABC-123"],
              "summary": "covalently targets BTK; AML and NSCLC excluded"}
        assert targets_in(tr, idx, id2) == ["BTK", "EGFR"], targets_in(tr, idx, id2)   # NSCLC/AML/Phase filtered
        # modality filter: pure radioligand/antibody dropped; small-molecule kept; SM+antibody COMBO kept
        assert _smallmol_program({"drugs": ["177Lu-PSMA-617"], "title": "PSMA radioligand"}, ["EGFR"]) == []
        assert _smallmol_program({"drugs": ["ABC-123"], "title": "EGFR inhibitor"}, ["EGFR"]) == ["EGFR"]
        assert _smallmol_program({"drugs": ["orelabrutinib", "anti-CD20 antibody"], "title": "BTK + antibody"}, ["BTK"]) == ["BTK"]
        assert _COVALENT_CLASS.search("Study of XYZ in KRAS G12C Mutation Advanced Solid Tumors")   # covalent-by-class
        cov = {"EGFR": {"NEOPHORE"}, "BTK": {"MERCK"}}
        assert cross_link({"sponsor": "NeoPhore Ltd"}, ["EGFR"], cov)["level"] == "LIKELY"
        assert cross_link({"sponsor": "Pfizer Inc"}, ["EGFR"], cov)["level"] == "POSSIBLE"
        assert cross_link({"sponsor": "X"}, ["KRAS"], cov)["level"] == "—"
        # 2c fingerprint: a covalent-design protocol scores high; a generic dose-escalation protocol doesn't
        assert fingerprint("Target occupancy in PBMC at 24/48h; covalent adduct by mass spectrometry")["score"] >= FP_THRESHOLD
        assert fingerprint("Maximum tolerated dose; pharmacokinetics; tumor response by RECIST")["score"] < FP_THRESHOLD
        # full-text-2c verdict: mechanism fires; bare word fires unless negated by comparator/non-covalent
        assert _pub_verdict("osimertinib is an irreversible EGFR-TKI")
        assert _pub_verdict("the inhibitor covalently binds Cys797 of EGFR")
        assert not _pub_verdict("pirtobrutinib is a non-covalent BTK inhibitor, unlike covalent BTK inhibitors")
        assert _ctis_date("GR: 29/06/2026") == "2026-06-29" and _ctis_date("n/a") == ""   # (a) CTIS date parse
        ph = lambda s: bool(re.search(r"\bphase i\b", s, re.I))                             # phase filter is \b, not substring
        assert ph("Phase I and Phase II (Integrated)") and not ph("Phase III") and not ph("Phase II")
        assert _novelty_tier(False, False, 0) == "NOVEL"                                    # (b) novelty tiers
        assert _novelty_tier(True, False, 0) == "EMERGING" and _novelty_tier(True, True, 5) == "KNOWN"
        # cross-registry dedup: same (firm,target) via NCT + CTIS -> one move (strongest covalency, earliest)
        dd = _dedup_moves([{"sponsor": "Verastem Inc.", "target": "KRAS", "covalency": "POSSIBLE", "posted": "2026-06-15"},
                           {"sponsor": "Verastem Inc", "target": "KRAS", "covalency": "LIKELY", "posted": "2026-06-09"}])
        assert len(dd) == 1 and dd[0]["covalency"] == "LIKELY", dd
        print("clinical_radar self-check OK")
    else:
        run(web="--web" in sys.argv, pubs="--pubs" in sys.argv, conf="--conf" in sys.argv)
