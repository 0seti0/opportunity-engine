"""Backtest the discovery thesis: does an as-of-2021 'emerging' signal predict 2022-2026 clinical emergence,
and does 'under-crowded' (Pharos TDL) sharpen it toward NOVEL targets?

v0 uses the reliable keyless legs: Europe PMC literature trajectory (signal) + CT.gov new trials (outcome,
non-circular) + Pharos TDL (crowding). The PATENT leg -- the real differentiator -- plugs into emergence()
once a PatentsView key exists. Point-in-time clean: signal reads only <=2021, outcome only >=2022.
ponytail: literature-as-signal first; term-search counts (gene mention != trial-of) are a known noise ceiling."""
import csv
import json
import os
import random
import statistics
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

CACHE = Path(__file__).parent / "backtest_cache.json"
_cache = json.loads(CACHE.read_text()) if CACHE.exists() else {}


def _get(url, method="GET", body=None, headers=None):
    key = method + url + (json.dumps(body) if body else "")   # token-free cache key -> stable across runs
    if key in _cache:
        return _cache[key]
    h = {"User-Agent": "Mozilla/5.0", "Accept": "application/json"}
    h.update(headers or {})
    data = json.dumps(body).encode() if body else None
    if body:
        h["Content-Type"] = "application/json"
    for attempt in range(4):                       # backoff on the flaky/rate-limited hosts; never swallow
        try:
            r = urllib.request.urlopen(urllib.request.Request(url, data=data, headers=h, method=method), timeout=45)
            out = json.loads(r.read())
            _cache[key] = out
            return out
        except urllib.error.HTTPError as e:
            if e.code in (429, 500, 502, 503) and attempt < 3:
                time.sleep(2 * (attempt + 1)); continue
            raise
        except Exception:                          # socket timeout / URLError -> retry with backoff (don't crash a long sweep)
            if attempt < 3:
                time.sleep(2 * (attempt + 1)); continue
            raise


def epmc_year(gene, year):                         # literature hitCount for a gene in one PUB_YEAR
    return _get("https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + urllib.parse.urlencode(
        {"query": f'{gene} AND PUB_YEAR:{year}', "format": "json", "resultType": "idlist", "pageSize": 1})).get("hitCount", 0)


def ct_trials(gene, start, end):                   # interventional trials mentioning gene, started in [start,end]
    return _get("https://clinicaltrials.gov/api/v2/studies?" + urllib.parse.urlencode(
        {"query.term": gene, "filter.advanced": f"AREA[StudyType]INTERVENTIONAL AND AREA[StartDate]RANGE[{start},{end}]",
         "countTotal": "true", "pageSize": 1, "format": "json"})).get("totalCount", 0)


def tdl(gene):                                     # Pharos Target Development Level (crowding); Tclin = has approved drug
    d = _get("https://pharos-api.ncats.io/graphql", "POST", {"query": f'{{target(q:{{sym:"{gene}"}}){{tdl}}}}'})
    return ((d.get("data") or {}).get("target") or {}).get("tdl", "")


LENS = "https://api.lens.org/patent/search"


def _lens_token():                                 # from LENS_TOKEN env or local .lens_token file (gitignored; never hardcode)
    return os.environ.get("LENS_TOKEN") or (
        (Path(__file__).parent / ".lens_token").read_text().strip()
        if (Path(__file__).parent / ".lens_token").exists() else None)


def lens_count(gene, start, end):                  # patents whose EARLIEST PRIORITY date in [start,end] mention gene
    tok = _lens_token()
    if not tok:
        raise RuntimeError("no Lens token: set LENS_TOKEN env or write layer1/.lens_token")
    body = {"query": {"bool": {
        "must": [{"query_string": {"query": gene, "fields": ["title", "abstract", "claim"]}}],
        "filter": [{"range": {"earliest_priority_claim_date": {"gte": start, "lte": end}}}]}}, "size": 0}
    return _get(LENS, "POST", body, headers={"Authorization": f"Bearer {tok}"}).get("total", 0)


def patent_emergence(gene):                        # SIGNAL (the real thesis): patent-filing growth into 2021, intent-weighted
    recent = lens_count(gene, "2019-01-01", "2021-12-31")
    base = lens_count(gene, "2015-01-01", "2018-12-31")
    return recent, (recent + 1) / (base + 1)       # ponytail: title/abstract/claim text-match; short symbols are noisy (the precision ceiling)


REPORTER = "https://api.reporter.nih.gov/v2/projects/search"


def reporter_year(gene, year):                     # NIH grants mentioning gene in ONE fiscal year (single-year totals reliable;
    body = {"criteria": {"advanced_text_search": {"operator": "and",                       # multi-year-list totals glitch -> sum singles)
            "search_field": "projecttitle,terms,abstracttext", "search_text": gene}, "fiscal_years": [year]},
            "include_fields": ["ProjectNum"], "limit": 1}
    return _get(REPORTER, "POST", body).get("meta", {}).get("total", 0)


def grant_emergence(gene):                         # SIGNAL (keyless): NIH-grant funding growth into the 2021 freeze
    recent = sum(reporter_year(gene, y) for y in (2019, 2020, 2021))
    base = sum(reporter_year(gene, y) for y in (2015, 2016, 2017, 2018))
    return recent, (recent + 1) / (base + 1)


def emergence(gene):                               # SIGNAL: literature growth INTO the 2021 freeze, only <=2021 data
    recent = epmc_year(gene, 2020) + epmc_year(gene, 2021)
    base = epmc_year(gene, 2016) + epmc_year(gene, 2017)
    return recent, (recent + 1) / (base + 1)       # ratio > 1 = rising into the freeze


def outcome(gene):                                 # OUTCOME: post-freeze clinical emergence = NEW trials 2022-2026
    return ct_trials(gene, "2022-01-01", "2026-12-31")


def universe(n, seed=7):
    rows = [r["symbol"] for r in csv.DictReader(open("/tmp/hgnc.tsv"), delimiter="\t")
            if r.get("locus_group") == "protein-coding gene"] or \
           [r["symbol"] for r in csv.DictReader(open("/tmp/hgnc.tsv"), delimiter="\t")]
    random.seed(seed)
    return random.sample(rows, n)


def run(n=40):
    genes = universe(n)
    rows = []
    for i, g in enumerate(genes):
        recent, emrg = emergence(g)
        rows.append({"gene": g, "recent_lit": recent, "emergence": emrg, "tdl": tdl(g), "new_trials": outcome(g)})
        if i % 10 == 9:
            CACHE.write_text(json.dumps(_cache))    # checkpoint the cache so a flaky run is resumable
    CACHE.write_text(json.dumps(_cache))

    cand = sorted((r for r in rows if r["recent_lit"] >= 5), key=lambda r: r["emergence"], reverse=True)  # min signal mass
    k = max(1, len(cand) // 3)
    top, bot = cand[:k], cand[-k:]
    mt = statistics.mean(r["new_trials"] for r in top)
    mb = statistics.mean(r["new_trials"] for r in bot)
    print(f"universe={n}  candidates(recent_lit>=5)={len(cand)}")
    print(f"mean NEW interventional trials 2022-26:  top-third emergence={mt:.1f}   bottom-third={mb:.1f}   "
          f"lift={mt / (mb + 1e-9):.1f}x")

    novel = [r for r in cand if r["tdl"] and r["tdl"] != "Tclin"]
    print("\nEMERGING x UNDER-CROWDED (high 2021 emergence, no approved drug) -> their REAL 2022-26 trial outcome:")
    for r in novel[:12]:
        hit = "EMERGED" if r["new_trials"] >= 5 else ("" if r["new_trials"] else "quiet")
        print(f"  {r['gene']:9} emergence={r['emergence']:4.1f}  recent_lit={r['recent_lit']:5}  {r['tdl']:6}"
              f" -> new_trials={r['new_trials']:4}  {hit}")


def run_patent(n=80):
    genes = universe(n)
    rows = []
    for i, g in enumerate(genes):
        prec, pemrg = patent_emergence(g)
        rows.append({"gene": g, "pat_recent": prec, "pat_emergence": pemrg, "tdl": tdl(g), "new_trials": outcome(g)})
        if i % 10 == 9:
            CACHE.write_text(json.dumps(_cache))
    CACHE.write_text(json.dumps(_cache))

    cand = sorted((r for r in rows if r["pat_recent"] >= 3), key=lambda r: r["pat_emergence"], reverse=True)  # patent activity = intent
    k = max(1, len(cand) // 3)
    top, bot = cand[:k], cand[-k:]
    mt = statistics.mean(r["new_trials"] for r in top)
    mb = statistics.mean(r["new_trials"] for r in bot)
    print(f"universe={n}  patent-active candidates(pat_recent>=3)={len(cand)}")
    print(f"mean NEW interventional trials 2022-26:  top-third PATENT-emergence={mt:.1f}   bottom-third={mb:.1f}   "
          f"lift={mt / (mb + 1e-9):.1f}x   (literature signal scored 0.0x on the same harness)")

    novel = [r for r in cand if r["tdl"] and r["tdl"] != "Tclin"]
    print("\nEMERGING(patents) x UNDER-CROWDED (no approved drug) -> their REAL 2022-26 trial outcome:")
    for r in novel[:12]:
        hit = "EMERGED" if r["new_trials"] >= 5 else ("" if r["new_trials"] else "quiet")
        print(f"  {r['gene']:9} pat_emrg={r['pat_emergence']:4.1f}  pat_recent={r['pat_recent']:4}  {r['tdl']:6}"
              f" -> new_trials={r['new_trials']:4}  {hit}")


def run_grant(n=80):
    genes = universe(n)
    rows = []
    for i, g in enumerate(genes):
        grec, gemrg = grant_emergence(g)
        rows.append({"gene": g, "grant_recent": grec, "grant_emergence": gemrg, "tdl": tdl(g), "new_trials": outcome(g)})
        if i % 10 == 9:
            CACHE.write_text(json.dumps(_cache))
    CACHE.write_text(json.dumps(_cache))

    cand = sorted((r for r in rows if r["grant_recent"] >= 3), key=lambda r: r["grant_emergence"], reverse=True)  # grant activity = intent
    k = max(1, len(cand) // 3)
    top, bot = cand[:k], cand[-k:]
    mt = statistics.mean(r["new_trials"] for r in top)
    mb = statistics.mean(r["new_trials"] for r in bot)
    print(f"universe={n}  grant-active candidates(grant_recent>=3)={len(cand)}")
    print(f"mean NEW interventional trials 2022-26:  top-third GRANT-emergence={mt:.1f}   bottom-third={mb:.1f}   "
          f"lift={mt / (mb + 1e-9):.1f}x   (literature signal = 0.0x on the same harness)")

    novel = [r for r in cand if r["tdl"] and r["tdl"] != "Tclin"]
    print("\nEMERGING(grants) x UNDER-CROWDED (no approved drug) -> their REAL 2022-26 trial outcome:")
    for r in novel[:12]:
        hit = "EMERGED" if r["new_trials"] >= 5 else ("" if r["new_trials"] else "quiet")
        print(f"  {r['gene']:9} grant_emrg={r['grant_emergence']:4.1f}  grant_recent={r['grant_recent']:4}  {r['tdl']:6}"
              f" -> new_trials={r['new_trials']:4}  {hit}")


def _spearman(xs, ys):                              # rank-correlation, stdlib (no scipy)
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        for pos, i in enumerate(order):
            r[i] = pos
        return r
    rx, ry, n = rank(xs), rank(ys), len(xs)
    mx, my = sum(rx) / n, sum(ry) / n
    cov = sum((rx[i] - mx) * (ry[i] - my) for i in range(n))
    vx = sum((rx[i] - mx) ** 2 for i in range(n)) ** 0.5
    vy = sum((ry[i] - my) ** 2 for i in range(n)) ** 0.5
    return cov / (vx * vy) if vx and vy else 0.0


def load_patent_counts(path):                       # BigQuery (gene,yr,n) CSV -> {gene: {yr: count}}
    d = {}
    for row in csv.DictReader(open(path)):
        d.setdefault(row["gene"], {})[int(row["yr"])] = int(row["n"])
    return d


def run_patent_bq():                                # the real patents-first backtest, on the cached BigQuery scan
    base_dir = Path(__file__).parent
    univ = json.loads((base_dir / "universe.json").read_text())
    pc = load_patent_counts(base_dir / "patent_counts.csv")
    rows = []
    for i, g in enumerate(univ):
        yrs = pc.get(g, {})
        p_base = sum(yrs.get(y, 0) for y in range(2015, 2019))
        p_recent = sum(yrs.get(y, 0) for y in range(2019, 2022))
        if p_recent < 10:                           # intent floor: gene must be patent-active by 2021 (>=10 keeps the network sweep tractable)
            continue
        rows.append({"gene": g, "p_base": p_base, "p_recent": p_recent,
                     "emergence": (p_recent + 1) / (p_base + 1), "tdl": tdl(g), "new_trials": outcome(g)})
        if i % 25 == 24:
            CACHE.write_text(json.dumps(_cache))
    CACHE.write_text(json.dumps(_cache))

    rows.sort(key=lambda r: r["emergence"], reverse=True)
    k = max(1, len(rows) // 3)
    top, bot = rows[:k], rows[-k:]
    mt = statistics.mean(r["new_trials"] for r in top)
    mb = statistics.mean(r["new_trials"] for r in bot)
    sp = _spearman([r["emergence"] for r in rows], [r["new_trials"] for r in rows])
    base_rate = statistics.mean(1 if r["new_trials"] >= 5 else 0 for r in rows)
    prec = statistics.mean(1 if r["new_trials"] >= 5 else 0 for r in rows[:20])
    print(f"patent-active candidates (p_recent>=10): {len(rows)}  (universe {len(univ)})")
    print(f"LIFT: top-third emergence new_trials={mt:.1f}  bottom-third={mb:.1f}  lift={mt / (mb + 1e-9):.1f}x"
          f"   (grants=1.3x, literature=0.0x)")
    print(f"  baseline check (median p_recent): top={statistics.median(r['p_recent'] for r in top):.0f}"
          f"  bottom={statistics.median(r['p_recent'] for r in bot):.0f}  (similar => lift isn't just an activity artifact)")
    print(f"Spearman(emergence, new_trials) = {sp:+.2f}")
    print(f"precision@20: {prec:.0%} of top-emergence genes EMERGED (>=5 new trials)  vs base rate {base_rate:.0%}")
    print("\nEMERGING(patents) x UNDER-CROWDED (not Tclin) -> their REAL 2022-26 trial outcome:")
    for r in [r for r in rows if r["tdl"] and r["tdl"] != "Tclin"][:15]:
        hit = "EMERGED" if r["new_trials"] >= 5 else ("" if r["new_trials"] else "quiet")
        print(f"  {r['gene']:9} emrg={r['emergence']:5.1f}  p_recent={r['p_recent']:4}  {r['tdl']:6}"
              f" -> new_trials={r['new_trials']:4}  {hit}")


if __name__ == "__main__":
    import sys
    mode = sys.argv[1] if len(sys.argv) > 1 else ("bq" if (Path(__file__).parent / "patent_counts.csv").exists()
                                                  else "patent" if _lens_token() else "grant")
    {"lit": run, "patent": run_patent, "grant": run_grant, "bq": run_patent_bq}[mode]()   # python3 backtest.py [bq|grant|patent|lit]
