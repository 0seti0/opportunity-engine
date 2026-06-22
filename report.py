"""Enriched opportunity report: each best-opportunity target with its EVIDENCE + PROVENANCE —
covalent-handle (CysDB/KLIFS), validation (CysDB drug-target flags + ClinicalTrials.gov), and any
covalent patents in our SureChEMBL feed (with attribution confidence). Every field is source-labelled.
"""
import csv
import json
import re
from collections import defaultdict
from pathlib import Path

import openpyxl

import backtest as bt
import window as W

HERE = Path(__file__).parent


def run():
    opps = json.load(open(HERE / "covalent_opportunities_final.json"))["best_opportunities"]
    feed = {e[0]: e for e in json.load(open(HERE / "covalent_warhead_feed_clean.json"))}
    tg = json.load(open(HERE / "patent_targets.json"))

    # covalent-feed patents attributed to each target (any confidence) -> the patent evidence/provenance
    pat = defaultdict(list)
    for pn, d in tg.items():
        e = feed.get(pn)
        if e and d.get("sym"):
            pat[d["sym"]].append({"patent": pn, "date": e[1], "assignee": W.norm_assignee(e[2].split(';')[0]),
                                  "title": e[3], "warheads": e[5], "attribution": d["conf"]})

    # UniProt + CysDB Fig4A-C drug-target flags
    uni = {}
    for r in csv.DictReader(open("/tmp/hgnc.txt"), delimiter="\t"):
        if (r.get("uniprot_ids") or "").strip():
            uni[r["symbol"]] = r["uniprot_ids"].split("|")[0].strip()
    cysflag = {}
    wb = openpyxl.load_workbook(str(HERE / "NIHMS1893018-supplement-2.xlsx"), read_only=True)
    it = wb["Fig4A-C"].iter_rows(values_only=True); next(it)
    for pid, fda, chembl, drugbank, lig in it:
        cysflag[pid] = (fda == "yes", chembl == "yes", drugbank == "yes", lig == "yes")

    out = []
    for o in opps:
        t = o["target"]
        f = cysflag.get(uni.get(t), (False, False, False, False))
        ey = W.earliest_clinical(t)                                  # ClinicalTrials.gov earliest interventional year
        total = bt.ct_trials(t, "2000-01-01", "2026-12-31")
        recent = bt.ct_trials(t, "2022-01-01", "2026-12-31")
        ps = sorted(pat.get(t, []), key=lambda x: x["date"], reverse=True)
        # only a TITLE/CLAIMS attribution is credible evidence the patent concerns this target; a lone
        # body pathway-mention is NOT (same bar as the covalent-closed exclusion). Keep discarded ones
        # in a labelled audit field so nothing is silently dropped.
        credible = [p for p in ps if p["attribution"] in ("title", "claims")]
        discarded = [p for p in ps if p["attribution"] == "body"]
        out.append({
            "rank": o["rank"], "target": t, "lane": o["lane"],
            "covalent_handle": o["handle"],
            "validation": {
                "fda_approved_drug_target": f[0], "chembl_target": f[1], "drugbank_target": f[2],
                "clinical_trials_ctgov": {"total_interventional": total, "recent_2022_26": recent,
                                          "earliest_entry_year": ey, "search_term": t},
            },
            "covalent_patents_surechembl": credible or "none — no title/claims covalent patent in feed (lane open)",
            "discarded_body_mentions": discarded or None,   # attributed by body pathway-mention only; NOT evidence
            "thesis": o["thesis"], "top_risk": o["top_risk"],
            "sources": ["SureChEMBL bulk 2026-06-15 (patents + biomedical_entities) — patents & attribution",
                        "ClinicalTrials.gov v2 API — trials", "CysDB / Cell Chem Biol 2023 — drug-target flags + scout-fragment cysteines",
                        "KLIFS — kinase pocket cysteines", "HGNC — gene/UniProt IDs",
                        "web verification (cited in covalent_bic_final.json) — existing covalent programs"],
        })
    bt.CACHE.write_text(json.dumps(bt._cache))
    json.dump(out, open(HERE / "opportunity_report.json", "w"), indent=1)

    for o in out:
        v = o["validation"]; ct = v["clinical_trials_ctgov"]
        flags = "".join(c for c, b in [("FDA", v["fda_approved_drug_target"]), ("ChEMBL", v["chembl_target"]),
                                        ("DrugBank", v["drugbank_target"])] if b) or "—"
        print(f"\n#{o['rank']} {o['target']}  [{o['lane']}]  handle: {o['covalent_handle']}")
        print(f"   validation: {flags} | CT.gov: {ct['total_interventional']} trials total, "
              f"{ct['recent_2022_26']} since 2022, earliest ~{ct['earliest_entry_year']}")
        pats = o["covalent_patents_surechembl"]
        if isinstance(pats, list):
            print(f"   covalent patents (title/claims) in feed: {len(pats)}")
            for p in pats[:3]:
                print(f"      - {p['patent']} ({p['date']}) {p['assignee']} | {p['title'][:46]} | {p['warheads']} | attr={p['attribution']}")
        else:
            print(f"   covalent patents in feed: {pats}")
        if o["discarded_body_mentions"]:
            d0 = o["discarded_body_mentions"][0]
            print(f"   [discarded — body pathway-mention, NOT evidence: {d0['patent']} {d0['assignee']} '{d0['title'][:34]}']")
    print("\nsaved opportunity_report.json")


if __name__ == "__main__":
    run()
