"""Structural-cysteine complement to the CysDB layer — fixes CysDB's false-negative on deep ATP-pocket
cysteines (RAF1 C424 hinge cysteine scored 'not ligandable' by scout-fragment chemoproteomics, yet it is
the canonical covalent-RAF handle = BRAF C532).

Source: KLIFS assigns every human kinase's 85 ATP-pocket residues a canonical position. Calibrated against
known covalent drugs: pos 52 = front-pocket 'roof' Cys (EGFR C797, BTK C481, JAK3 C909, HER2 C805);
pos 46-48 = hinge Cys (BRAF C532, RAF1 C424); pos 4-9 = P-loop Cys (YES1 C287, FGFR1 C488); pos 45 = gatekeeper.

Opportunity = kinase with a druggable-position pocket cysteine, FDA-validated, covalent-OPEN (no covalent
patent in our feed), non-essential. Flags which were MISSED by CysDB scout-fragments (the recovery set).
"""
import csv
import json
from pathlib import Path

import openpyxl

HERE = Path(__file__).parent
# KLIFS pocket position -> covalent-handle class (1-based in the 85-residue pocket)
REGION = {}
for i in (4, 5, 6, 7, 8, 9): REGION[i] = "P-loop"
REGION[45] = "gatekeeper"
for i in (46, 47, 48): REGION[i] = "hinge"
for i in (49, 50, 51, 52): REGION[i] = "front-pocket"   # 52 = the EGFR-C797/BTK-C481 roof position
RANK = {"front-pocket": 0, "hinge": 1, "P-loop": 2, "gatekeeper": 3}


def run():
    klifs = json.load(open("/tmp/klifs.json"))
    # CysDB Fig4A-C: uniprot -> drug-target + scout-fragment-ligandable flags
    cys = {}
    wb = openpyxl.load_workbook(str(HERE / "NIHMS1893018-supplement-2.xlsx"), read_only=True)
    it = wb["Fig4A-C"].iter_rows(values_only=True); next(it)
    for pid, fda, chembl, drugbank, lig in it:
        if pid:
            cys[pid] = {"fda": fda == "yes", "chembl": chembl == "yes",
                        "drugbank": drugbank == "yes", "cysdb_ligandable": lig == "yes"}
    # covalent-patent set (title/claims only — a body-mention is too weak to mark "covalent-closed",
    # the RAF1 false-positive lesson) + pan-essential set (same as cysdb.py)
    tg = json.load(open(HERE / "patent_targets.json"))
    covset = {d["sym"] for d in tg.values() if d.get("sym") and d.get("conf") in ("title", "claims")}
    cov_body_only = {d["sym"] for d in tg.values() if d.get("sym")} - covset
    ess = {ln.split("\t")[0] for ln in open("/tmp/ceg2.txt").read().splitlines()[1:] if ln.strip()}
    ess |= {ln.split(" (")[0].strip() for ln in open("/tmp/depmap_common_ess.csv").read().splitlines()[1:] if ln.strip()}

    rows = []
    for k in klifs:
        g, pocket = k.get("HGNC"), k.get("pocket") or ""
        if not g or len(pocket) != 85:
            continue
        handles = sorted({(REGION[i + 1], i + 1) for i, a in enumerate(pocket) if a == "C" and (i + 1) in REGION},
                         key=lambda x: (RANK[x[0]], x[1]))
        if not handles:
            continue
        c = cys.get(k.get("uniprot"), {})
        rows.append({"sym": g, "uniprot": k.get("uniprot"), "group": k.get("group"),
                     "handles": [f"{cl}@{p}" for cl, p in handles], "best": handles[0][0],
                     "fda": c.get("fda", False), "chembl": c.get("chembl", False),
                     "drugbank": c.get("drugbank", False), "cysdb_ligandable": c.get("cysdb_ligandable", False),
                     "covalent_patented": g in covset, "essential": g in ess})

    # dedupe by symbol (KLIFS lists isoforms/duplicates)
    bysym = {}
    for r in rows:
        b = bysym.setdefault(r["sym"], r)
        if RANK[r["best"]] < RANK[b["best"]]:
            bysym[r["sym"]] = r
    rows = list(bysym.values())

    opps = [r for r in rows if r["fda"] and not r["covalent_patented"] and not r["essential"]]
    opps.sort(key=lambda r: (RANK[r["best"]], not r["cysdb_ligandable"]))
    recovered = [r for r in opps if not r["cysdb_ligandable"]]   # structure found a handle CysDB missed
    json.dump({"opportunities": opps, "all_pocket_cys_kinases": len(rows)},
              open(HERE / "structural_cys_opportunities.json", "w"), indent=1)

    print(f"human kinases with a druggable-position pocket cysteine: {len(rows)}")
    resurfaced = sorted(r["sym"] for r in opps if r["sym"] in cov_body_only)
    print(f"  FDA-validated, covalent-open, non-essential: {len(opps)}")
    print(f"  ...of which CysDB scout-fragments MISSED (recovered by structure): {len(recovered)}")
    print(f"  ...re-surfaced by the title/claims fix (were body-only false-closed): {len(resurfaced)} {resurfaced}")
    print(f"  RAF1 now an opportunity? {'YES' if any(r['sym']=='RAF1' for r in opps) else 'no'}"
          f"  (handles={next((r['handles'] for r in opps if r['sym']=='RAF1'), '—')})")
    # sanity: known covalent-drugged kinases must be EXCLUDED (covalent_patented)
    print("\nsanity — known covalent-drugged kinases (should be covalent_patented=excluded):")
    for g in ("EGFR", "BTK", "ERBB2", "JAK3"):
        r = bysym.get(g)
        if r: print(f"  {g:6} handles={r['handles']} covalent_patented={r['covalent_patented']}")
    print("\nTOP structural-cysteine opportunities (front-pocket & hinge first; ✚ = CysDB missed it):")
    for r in opps[:30]:
        tags = "".join(t for t, f in [("F", r["fda"]), ("C", r["chembl"]), ("D", r["drugbank"])] if f)
        rec = "✚" if not r["cysdb_ligandable"] else " "
        print(f"  {rec} {r['sym']:9} {','.join(r['handles']):24} [{tags}] {r['group']}")


if __name__ == "__main__":
    run()
