"""Walking skeleton: CT.gov trials + openFDA approvals -> extract gene/drug mentions -> resolve (both resolvers)
-> append-only claim log. Two evidence sources, both paths. Reuses identity + drug."""
import datetime
import hashlib
import json
import re
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

from identity import Outcome, build_index, load_rows, resolve
from drug import resolve_drug

LOG = Path(__file__).parent / "claims.jsonl"


def _get(url):
    return json.load(urllib.request.urlopen(
        urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"}), timeout=60))


def fetch_trials(condition, n):
    return _get("https://clinicaltrials.gov/api/v2/studies?" + urllib.parse.urlencode({
        "query.cond": condition, "filter.advanced": "AREA[StudyType]INTERVENTIONAL",
        "pageSize": n, "format": "json"})).get("studies", [])


def fetch_labels(indication, n):
    return _get("https://api.fda.gov/drug/label.json?" + urllib.parse.urlencode({
        "search": f'indications_and_usage:"{indication}"', "limit": n})).get("results", [])


def fetch_patents(query, n):               # Google Patents XHR (keyless); priority_date = earliest pre-trial signal
    # ponytail: undocumented endpoint, one page; swap to PatentsView (needs free API key) if it breaks/rate-limits
    url = "https://patents.google.com/xhr/query?" + urllib.parse.urlencode(
        {"url": f"q={query}&type=PATENT&num={n}&sort=new"})   # newest filings -> recent priority dates surface
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    cl = json.load(urllib.request.urlopen(req, timeout=60)).get("results", {}).get("cluster", [])
    return [r["patent"] for r in (cl[0].get("result", []) if cl else [])]


SC = "https://www.surechembl.org/api"


def _sc_title(meta):                       # English title -> Path A text
    for t in meta.get("titles", []):
        if t.get("lang") == "en":
            return " ".join(t.get("titles", []))
    return ""


def fetch_patents_sc(query, n):            # SureChEMBL (EBI, keyless, global). No server date-sort -> pull a page, take newest
    url = SC + "/search/content?" + urllib.parse.urlencode({"query": query, "page": 1, "itemsPerPage": 50})
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0", "Accept": "application/json"}, method="POST")
    docs = json.load(urllib.request.urlopen(req, timeout=60))["data"]["results"]["documents"]
    docs.sort(key=lambda d: d.get("metadata", {}).get("pd", ""), reverse=True)
    return docs[:n]


EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest"
EPMC_ANN = "https://www.ebi.ac.uk/europepmc/annotations_api/annotationsByArticleIds"


def fetch_literature(query, n):            # Europe PMC (keyless), newest-first; papers are often the earliest signal
    return _get(EPMC + "/search?" + urllib.parse.urlencode({
        "query": query, "format": "json", "resultType": "core",
        "pageSize": n, "sort": "P_PDATE_D desc"})).get("resultList", {}).get("result", [])


def gene_mentions(article_id):             # text-mined gene SPANS (recognition) in title/abstract -> OUR HGNC index normalizes
    # ponytail: Europe PMC's own normalization tags are cross-species junk; we keep only the `exact` span and resolve it ourselves
    arts = _get(EPMC_ANN + "?" + urllib.parse.urlencode(
        {"articleIds": article_id, "type": "Gene_Proteins", "format": "JSON"}))
    anns = arts[0].get("annotations", []) if arts else []
    return {a.get("exact", "") for a in anns if a.get("section", "").startswith(("Title", "Abstract"))}


_appdate = {}


def _ymd(s):                               # normalize YYYYMMDD -> YYYY-MM-DD; leave already-hyphenated dates
    return f"{s[:4]}-{s[4:6]}-{s[6:8]}" if len(s) == 8 and s.isdigit() else s


def approval_date(app_no):                 # real original-approval date from Drugs@FDA, normalized; "" if unknown
    if not app_no:
        return ""
    if app_no not in _appdate:
        try:
            res = _get("https://api.fda.gov/drug/drugsfda.json?" + urllib.parse.urlencode(
                {"search": f'application_number:"{app_no}"', "limit": 1})).get("results", [])
            ap = [s["submission_status_date"] for s in (res[0].get("submissions", []) if res else [])
                  if s.get("submission_status") == "AP" and s.get("submission_status_date")]
            _appdate[app_no] = _ymd(min(ap)) if ap else ""      # earliest AP submission = first approval
        except Exception:                  # drugsfda miss/404 -> unknown; claim still records as approval evidence
            _appdate[app_no] = ""
    return _appdate[app_no]


# ponytail: function-word stoplist kills "FOR"->WWOX-alias noise from all-caps patent titles. NOT a stoplist of real genes
# that happen to be words (SET/MET/REST stay) -- the real fix is a Path A precision/NER layer, this is the cheap guard.
STOP = frozenset("FOR AND THE NOT ALL ANY MAY CAN HAS HAD WAS ARE BUT NEW USE VIA PER OUR OFF OUT WHO HOW NOW OWN WAY ITS ONE TWO NON PRE".split())


def genes_in(text):                        # uppercase symbol-like tokens; the resolver rejects non-genes
    return {t for t in re.findall(r"\b[A-Z][A-Z0-9]{1,6}\b", text) if t not in STOP}


def inclusion_only(eligibility):           # keep only the inclusion half -> drops "no prior anti-TNF" exclusion false positives
    return re.split(r"(?i)exclusion criteria", eligibility)[0]


def make_claim(subject, source, evidence_type, source_id, valid_time, actor, path):
    # ponytail: claim_key omits content_hash -> no revised-record versioning yet
    key = hashlib.sha1(f"{source}|{source_id}|{evidence_type}|{subject}".encode()).hexdigest()[:16]
    return {"subject": subject, "evidence_type": evidence_type, "source": source,
            "source_record_id": source_id, "actor": actor, "valid_time": valid_time,
            "observed_at": str(datetime.date.today()), "path": path, "claim_key": key}


def run():
    gidx = build_index(load_rows())
    seen = {json.loads(l)["claim_key"] for l in LOG.read_text().splitlines()} if LOG.exists() else set()
    drug_cache, claims = {}, []

    def hits(text, drugs):                 # text -> genes (Path A) ; drugs -> targets (Path B)
        out = []
        for m in genes_in(text):
            r = resolve(m, gidx)
            if r.outcome is Outcome.RESOLVED:
                out.append((r.hgnc_id, "A"))
        for d in drugs:
            if d not in drug_cache:
                drug_cache[d] = resolve_drug(d, gidx)
            out += [(hid, "B") for hid in drug_cache[d].targets]
        return out

    for s in fetch_trials("HER2-positive breast cancer", 5):       # evidence_type: trial
        p = s["protocolSection"]
        text = (p["identificationModule"].get("briefTitle", "") + " " +
                " ; ".join(p.get("conditionsModule", {}).get("conditions", []) or []) + " " +
                inclusion_only(p.get("eligibilityModule", {}).get("eligibilityCriteria", "")))
        drugs = [i["name"] for i in p.get("armsInterventionsModule", {}).get("interventions", []) or []
                 if i.get("type") in ("DRUG", "BIOLOGICAL")]
        nct = p["identificationModule"]["nctId"]
        vt = _ymd(p.get("statusModule", {}).get("studyFirstPostDateStruct", {}).get("date", ""))
        actor = p.get("sponsorCollaboratorsModule", {}).get("leadSponsor", {}).get("name", "")
        for hid, path in hits(text, drugs):
            claims.append(make_claim(hid, "ctgov", "trial", nct, vt, actor, path))

    for lab in fetch_labels("breast cancer", 10):                  # evidence_type: approval
        of = lab.get("openfda", {})
        drug = (of.get("generic_name") or [""])[0]
        vt = approval_date((of.get("application_number") or [""])[0])    # real approval date, not label effective_time
        text = " ".join(lab.get("indications_and_usage", []) or []) + " " + " ".join(lab.get("mechanism_of_action", []) or [])
        for hid, path in hits(text, [drug] if drug else []):
            claims.append(make_claim(hid, "openfda", "approval", lab.get("set_id", ""),
                                     vt, (of.get("manufacturer_name") or [""])[0], path))

    for pat in fetch_patents("breast cancer", 12):                 # evidence_type: patent (earliest, pre-trial)
        text = re.sub(r"<[^>]+>", " ", pat.get("title", "") + " " + pat.get("snippet", ""))
        vt = _ymd(pat.get("priority_date") or pat.get("filing_date") or "")
        for hid, path in hits(text, []):                           # Path A only; patents have no clean structured drug field
            claims.append(make_claim(hid, "patents", "patent", pat.get("publication_number", ""),
                                     vt, pat.get("assignee", ""), path))

    for pat in fetch_patents_sc("breast cancer", 12):              # 2nd patent source: SureChEMBL (global, EBI). Same "patent" type
        meta = pat.get("metadata", {})
        vt = _ymd(meta.get("pd", ""))
        for hid, path in hits(_sc_title(meta), []):                # Path A (title only); chemistry->target Path B is a later upgrade
            claims.append(make_claim(hid, "surechembl", "patent", meta.get("pn") or pat.get("docId", ""),
                                     vt, pat.get("pa", ""), path))

    for art in fetch_literature("breast cancer therapeutic target", 15):   # evidence_type: literature (earliest signal)
        aid = f"{art.get('source')}:{art.get('id')}"
        vt = art.get("firstPublicationDate", "")                # already YYYY-MM-DD
        actor = art.get("journalInfo", {}).get("journal", {}).get("title", "") or art.get("source", "")
        for m in gene_mentions(aid):
            r = resolve(m, gidx)                                # NER recognition + our human-HGNC normalization (3-outcome)
            if r.outcome is Outcome.RESOLVED:
                claims.append(make_claim(r.hgnc_id, "europepmc", "literature", aid, vt, actor, "A"))

    new = 0
    with open(LOG, "a") as log:
        for c in claims:
            if c["claim_key"] not in seen:
                seen.add(c["claim_key"]); log.write(json.dumps(c) + "\n"); new += 1
    print(f"{new} new claims -> {LOG.name} ({len(seen)} total)")

    types = defaultdict(set)               # the payoff: targets backed by >1 independent evidence type
    for l in LOG.read_text().splitlines():
        c = json.loads(l)
        types[c["subject"]].add(c["evidence_type"])
    for hid, ts in sorted(types.items()):
        if len(ts) > 1:
            print(f"  CORROBORATED {hid}: {sorted(ts)}")


if __name__ == "__main__":
    run()
