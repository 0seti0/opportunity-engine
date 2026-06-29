"""Leg 2 — clinical early-warning radar. Catches covalent moves BEFORE the patent publishes, from
ClinicalTrials.gov (no chemistry). New Phase-1/Early-Phase-1 industry trials -> classify covalency by
  (Step 1)  diff new trials
  (Step 2a) patent cross-link to our covalent census (incl. the Day-1 radar top-up)
  (Step 2b) live web check on the drug code-name (claude -p, no API key; optional)
The clinical-DESIGN fingerprint (2c) and the unify/cross-link to Leg 1 are the next pieces.

  python clinical_radar.py            # diff new Ph1 trials + patent cross-link (fast)
  python clinical_radar.py --web      # also run the code-name web check (slow: claude -p per drug)
  python clinical_radar.py --selftest # offline checks
"""
import json
import re
import sys
import urllib.parse
import urllib.request
from collections import defaultdict
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


def _smallmol_program(tr, targets):
    """Keep only covalently-tractable targets of a small-molecule (non-radioligand/CAR-T/antibody/probe) trial."""
    dt = " ".join(tr["drugs"] + [tr["title"]])
    if _NONSM.search(dt) or _PROBE.search(dt):
        return []
    return [t for t in targets if covalently_tractable(t)]


def new_trials(since=None, page_size=200):
    """Newest-first industry Phase-1/Early-Phase-1 interventional trials (stop at `since` if given)."""
    q = {"filter.advanced": "AREA[Phase](PHASE1 OR EARLY_PHASE1) AND AREA[StudyType]INTERVENTIONAL "
                            "AND AREA[LeadSponsorClass]INDUSTRY",
         "sort": "StudyFirstPostDate:desc", "pageSize": page_size, "format": "json"}
    out = []
    for s in bt._get(f"{CTGOV}?{urllib.parse.urlencode(q)}").get("studies", []):
        ps = s.get("protocolSection", {}); idm = ps.get("identificationModule", {})
        posted = ps.get("statusModule", {}).get("studyFirstPostDateStruct", {}).get("date", "")
        if since and posted and posted < since:
            break
        out.append({"nct": idm.get("nctId"), "posted": posted,
                    "sponsor": ps.get("sponsorCollaboratorsModule", {}).get("leadSponsor", {}).get("name", ""),
                    "title": idm.get("briefTitle", ""),
                    "drugs": [i.get("name") for i in ps.get("armsInterventionsModule", {}).get("interventions", [])
                              if i.get("type") == "DRUG" and "placebo" not in (i.get("name") or "").lower()],
                    "summary": (ps.get("descriptionModule", {}).get("briefSummary", "") or "")[:800],
                    "design_text": _design_text(ps)})
    return out


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
    text = " ".join([trial["title"]] + trial["drugs"] + [trial["summary"]])
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
              f"\"evidence\": \"\"}}  (covalent null if you can't tell).")
    try:
        r = subprocess.run(["claude", "-p", prompt, "--output-format", "json"],
                           capture_output=True, text=True, timeout=300)
        o = _extract_json(json.loads(r.stdout).get("result", ""))
        return o if isinstance(o, dict) else {}
    except Exception:
        return {}


def run(web=False, pubs=False):
    seen = set(json.loads(SEEN.read_text())) if SEEN.exists() else set()
    idx, id2sym = gene_index()
    covidx = covalent_index()
    trials = new_trials()
    events = []
    for tr in trials:
        if tr["nct"] in seen or not tr["drugs"]:
            continue                         # seen, or a diagnostic/imaging/device trial (no DRUG intervention)
        targets = _smallmol_program(tr, targets_in(tr, idx, id2sym))
        if not targets:
            continue                         # no covalently-tractable small-molecule target -> 2a/2b can't classify
        xl = cross_link(tr, targets, covidx)
        ev = {**tr, "targets": targets, "covalency": xl["level"], "why": xl["why"], "target": xl["target"]}
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
            if w.get("covalent") is True:
                ev["covalency"], ev["why"] = "LIKELY", f"web: {(w.get('evidence') or '')[:80]}"
        ev.pop("design_text", None)                          # drop the bulky text from the event
        events.append(ev)

    rank = {"LIKELY": 0, "FINGERPRINT": 1, "POSSIBLE": 2, "—": 3}
    events.sort(key=lambda e: rank.get(e["covalency"], 3))
    for e in events:
        if e["covalency"] != "—":
            print(f"  [{e['covalency']:8}] {e['nct']} {e['posted']} {firm(e['sponsor'])[:20]:20} "
                  f"{(e['drugs'][0] if e['drugs'] else '?')[:22]:22} {e['target']}  ({e['why']})")
    SEEN.write_text(json.dumps(sorted(seen | {t['nct'] for t in trials if t['nct']})))
    json.dump(events, open(HERE / "clinical_events.json", "w"), indent=1)
    n_cov = sum(1 for e in events if e["covalency"] != "—")
    print(f"\nclinical radar: {len(trials)} new Ph1 trials · {len(events)} with a resolvable target · "
          f"{n_cov} covalent-suspect (cross-link{'+web' if web else ''}).")
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
        # modality filter: a radioligand/antibody trial is dropped; a small-molecule inhibitor kept
        assert _smallmol_program({"drugs": ["177Lu-PSMA-617"], "title": "PSMA radioligand"}, ["EGFR"]) == []
        assert _smallmol_program({"drugs": ["ABC-123"], "title": "EGFR inhibitor"}, ["EGFR"]) == ["EGFR"]
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
        print("clinical_radar self-check OK")
    else:
        run(web="--web" in sys.argv, pubs="--pubs" in sys.argv)
