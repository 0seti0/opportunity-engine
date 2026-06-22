"""CLEAN, PRE-REGISTERED backtest v2 (panel-prescribed, freeze shifted EARLIER so the mechanism-resolved
outcome actually exists in curated data). Does PATENT-emergence predict a target's FIRST mechanism-targeting
clinical entry, among targets that were trial-naive at the freeze?

WHY v2: v1 froze at end-2021 and got a 0% base rate — Open Targets has not yet curated 2022-26 first-entries
(INHBE's 2025 FIH shows as None). Fix: freeze at 2017 so the outcome window 2018-2022 is FULLY curated
(verified: KRAS's first mechanism-trial = 2018-08-27 is present). This also lets us define establishment
POINT-IN-TIME from Open Targets trial dates instead of anachronistic current Pharos TDL.

================= PRE-REGISTRATION (frozen BEFORE looking at the outcome) =================
FREEZE: end-2017.
UNIVERSE: random draw of HGNC protein-coding genes (symbol len>=4, not in BLOCKLIST, has Ensembl id),
  EXCLUDING any gene with a mechanism-targeting interventional trial STARTED at/before 2017
  [Open Targets drug->target->trialStartDate, point-in-time] => trial-naive-at-freeze population only.
SIGNAL: patent-emergence = (patents 2015-2017 + 1) / (patents 2012-2014 + 1), de-noised abstract counts (BigQuery).
OUTCOME (binary): gene acquired its FIRST mechanism-targeting interventional trial in 2018..2022
  (zero at/before 2017 by construction). Resolved by Open Targets MoA, NOT free-text symbol.
ESTIMATOR: Spearman(signal, outcome) + precision@20 vs base rate, with a LABEL-PERMUTATION p (one-sided).
DECISION RULE: patent-emergence SURVIVES iff Spearman POSITIVE with permutation p < 0.05.
  If null/negative -> "patent-emergence is NOT a demonstrated predictor" (report it; do not retune).
KNOWN LIMITATION: even 2018-2022 has mild Open Targets curation lag at the tail; treat as mildly conservative.
GUARDRAILS: one universe def, one outcome, one estimator, fixed here. BLOCKLIST + len filter are pre-freeze,
  outcome-independent. If fragile across the sample, report fragile, not positive.
==========================================================================================
"""
import csv
import json
import random
import statistics
from pathlib import Path

import backtest as bt                                   # reuse _get, _cache, CACHE, load_patent_counts, _spearman

OT = "https://api.platform.opentargets.org/api/v4/graphql"
HERE = Path(__file__).parent
BLOCKLIST = {"SARS2", "CASK", "PTPRC", "ITGAE", "AMN1", "REST", "CAMP", "GANC", "MARS1", "LARS1", "IARS1"}

FREEZE = 2017
BASE_YEARS = range(2012, 2015)                          # 2012-2014
RECENT_YEARS = range(2015, 2018)                        # 2015-2017
OUTCOME_LO, OUTCOME_HI = 2018, 2022                     # fully curated in Open Targets


def ensembl_map():
    return {r["symbol"]: (r.get("ensembl_gene_id") or "")
            for r in csv.DictReader(open("/tmp/hgnc.tsv"), delimiter="\t")}


def ot_first_trial_year(ensembl):                       # earliest mechanism-targeting trial YEAR (Open Targets), or None
    q = ('{ target(ensemblId:"%s"){ drugAndClinicalCandidates{ rows{ clinicalReports{ trialStartDate } } } } }'
         % ensembl)
    d = bt._get(OT, "POST", {"query": q})
    rows = (((d.get("data") or {}).get("target") or {}).get("drugAndClinicalCandidates") or {}).get("rows") or []
    dates = [cr.get("trialStartDate") for r in rows for cr in (r.get("clinicalReports") or []) if cr.get("trialStartDate")]
    return min(int(x[:4]) for x in dates) if dates else None


def patent_emergence(g, pc):
    y = pc.get(g, {})
    base = sum(y.get(yr, 0) for yr in BASE_YEARS)
    recent = sum(y.get(yr, 0) for yr in RECENT_YEARS)
    return recent, (recent + 1) / (base + 1)


def _perm_p(signal, outcome, observed, k=2000):         # one-sided label-permutation p for Spearman>=observed
    rng = random.Random(11)
    ge = 0
    for _ in range(k):
        sh = outcome[:]
        rng.shuffle(sh)
        if bt._spearman(signal, sh) >= observed:
            ge += 1
    return (ge + 1) / (k + 1)


def run(target_n=400, screen_cap=800):
    emap = ensembl_map()
    pc = bt.load_patent_counts(HERE / "patent_counts.csv")
    pool = [g for g in json.loads((HERE / "universe.json").read_text())
            if g not in BLOCKLIST and emap.get(g)]
    random.Random(7).shuffle(pool)

    rows, screened, excl_pretrial = [], 0, 0
    for g in pool:
        if len(rows) >= target_n or screened >= screen_cap:
            break
        screened += 1
        fy = ot_first_trial_year(emap[g])
        if fy is not None and fy <= FREEZE:             # had a mechanism trial at/before freeze -> not trial-naive
            excl_pretrial += 1
            continue
        p_recent, p_emrg = patent_emergence(g, pc)
        rows.append({"gene": g, "p_recent": p_recent, "emergence": p_emrg,
                     "outcome": 1 if (fy is not None and OUTCOME_LO <= fy <= OUTCOME_HI) else 0, "first_year": fy})
        if screened % 25 == 0:
            bt.CACHE.write_text(json.dumps(bt._cache))
    bt.CACHE.write_text(json.dumps(bt._cache))

    sig = [r["emergence"] for r in rows]
    out = [r["outcome"] for r in rows]
    sp = bt._spearman(sig, out)
    p = _perm_p(sig, out, sp)
    base_rate = statistics.mean(out) if out else 0
    rows.sort(key=lambda r: r["emergence"], reverse=True)
    prec = statistics.mean(r["outcome"] for r in rows[:20])
    print(f"freeze={FREEZE}  screened={screened}  excluded not-trial-naive={excl_pretrial}")
    print(f"CLEAN universe (trial-naive at {FREEZE}): N={len(rows)}  |  emerged (first trial {OUTCOME_LO}-{OUTCOME_HI}): "
          f"{sum(out)} ({base_rate:.0%} base rate)")
    print(f"Spearman(patent-emergence, first-entry) = {sp:+.3f}   permutation p = {p:.3f}")
    print(f"precision@20 (top patent-emergence) = {prec:.0%}  vs base {base_rate:.0%}")
    verdict = "SURVIVES" if (sp > 0 and p < 0.05) else "FAILS (not a demonstrated predictor)"
    print(f"PRE-REGISTERED VERDICT: patent-emergence {verdict}")
    print("\ntop patent-emergence in the clean universe -> did they emerge 2018-2022?")
    for r in rows[:15]:
        print(f"  {r['gene']:9} emrg={r['emergence']:5.1f}  p_recent={r['p_recent']:4}"
              f"  first_trial={r['first_year']}  outcome={r['outcome']}")


if __name__ == "__main__":
    run()
