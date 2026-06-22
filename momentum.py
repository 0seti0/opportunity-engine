"""Protocol A — HOT EMERGING OPEN covalent targets. Over the covalently-tractable, covalent-OPEN,
non-essential, druggable universe (emerging_universe.json), score RECENT momentum from the newest data:
  accel  = annualized literature rate 2025-26 vs 2021-23 baseline (Europe PMC; includes preprints)
  preprints_25_26 = bioRxiv/medRxiv signal (earliest lead)   new_trials_24_26 = ClinicalTrials.gov
Momentum is a TIMING/prioritization signal over the structurally-open universe — NOT a prediction of
drug success (that was refuted). Surfaces which open covalent lanes are heating up NOW.
"""
import csv
import json
import time
import urllib.parse
from pathlib import Path

import openpyxl

import backtest as bt

HERE = Path(__file__).parent


def epmc(q):
    return bt._get("https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + urllib.parse.urlencode(
        {"query": q, "format": "json", "resultType": "idlist", "pageSize": 1})).get("hitCount", 0)


def rng(g, y0, y1):
    return epmc(f'{g} AND PUB_YEAR:[{y0} TO {y1}]')


def run():
    U = json.load(open(HERE / "emerging_universe.json"))["drug_universe"]
    uni2sym = {}
    for r in csv.DictReader(open("/tmp/hgnc.txt"), delimiter="\t"):
        for u in (r.get("uniprot_ids") or "").split("|"):
            if u.strip():
                uni2sym[u.strip()] = r["symbol"]
    fda = set()
    wb = openpyxl.load_workbook(str(HERE / "NIHMS1893018-supplement-2.xlsx"), read_only=True)
    it = wb["Fig4A-C"].iter_rows(values_only=True); next(it)
    for pid, f, c, d, l in it:
        if f == "yes" and uni2sym.get(pid):
            fda.add(uni2sym[pid])

    rows = []
    for i, g in enumerate(U):
        recent, base = rng(g, 2025, 2026), rng(g, 2021, 2023)
        rows.append({"sym": g, "recent": recent, "base": base,
                     "accel": round((recent / 2) / (base / 3 + 1), 2), "fda": g in fda})
        if i % 50 == 49:
            bt.CACHE.write_text(json.dumps(bt._cache)); print(f"  scanned {i+1}/{len(U)}", flush=True)
        time.sleep(0.05)
    bt.CACHE.write_text(json.dumps(bt._cache))

    hot = sorted((r for r in rows if r["recent"] >= 25), key=lambda r: -r["accel"])
    for r in hot[:40]:                                   # enrich the top with earliest + clinical confirm
        r["preprints_25_26"] = epmc(f'{r["sym"]} AND PUB_YEAR:[2025 TO 2026] AND (SRC:PPR)')
        r["new_trials_24_26"] = bt.ct_trials(r["sym"], "2024-01-01", "2026-12-31")
    bt.CACHE.write_text(json.dumps(bt._cache))
    json.dump({"all": rows, "hot_ranked": hot}, open(HERE / "emerging_open_ranked.json", "w"), indent=1)

    print("\n=== HOT EMERGING OPEN COVALENT TARGETS (covalent-open + ligandable, ranked by literature acceleration) ===")
    print(f"{'target':10}{'accel':>6}{'recent':>7}{'base':>6}{'preprint':>9}{'newtrials':>10}  status")
    for r in hot[:30]:
        print(f"{r['sym']:10}{r['accel']:>6}{r['recent']:>7}{r['base']:>6}"
              f"{r.get('preprints_25_26','-'):>9}{r.get('new_trials_24_26','-'):>10}  {'VALIDATED' if r['fda'] else 'emerging'}")
    print(f"\nscanned {len(rows)} | {len(hot)} with real recent activity | saved emerging_open_ranked.json")


if __name__ == "__main__":
    run()
