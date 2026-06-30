"""CovalentInDB 2.0 — deterministic covalent PRIOR-ART layer (static dataset, Jul-2024 cutoff).

Per HUMAN target gene: is there known covalent chemistry, which warhead(s), which bonded residue(s)
(e.g. EGFR CYS-797), the best kinact/Ki, and how many distinct chemotypes. This REPLACES the non-deterministic
`claude -p` guess for the *chemistry-precedent* question — but it is explicitly NOT a clinical-stage or recency
source: CovInDB_All.csv has no phase/approval field, and recent programs are absent (verified: WRN/VVD-214 = 0
rows, the canonical TRAP). So the dossier pairs this with the CT.gov clinical leg (stage) and the citation-
guarded `claude -p` (resistance / post-2024 gaps). Raw CSV is gitignored (IP-sensitive, downloaded once).
"""
import csv
import json
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).parent
CSV = HERE / "CovInDB_All.csv"
CACHE = HERE / ".covindb.json"


def _num(s):
    try:
        return float((s or "").replace(",", ""))
    except ValueError:
        return None


def index():
    """gene(UPPER) -> {n_inhibitors, n_chemotypes, warheads[], residues[], cys_sites[], best_kinact_ki}.
    Human-only (drops SARS-CoV-2/rat/parasite rows). Cached to .covindb.json (rebuild: delete the cache)."""
    if CACHE.exists():
        return json.load(open(CACHE))
    if not CSV.exists():
        return {}                                        # raw data not downloaded -> empty index; dossier shows "none"
    agg = defaultdict(lambda: {"n": 0, "chemo": set(), "warheads": set(), "residues": set(), "kk": None})
    for r in csv.DictReader(open(CSV)):
        if "Homo sapiens" not in r["Target_Taxonomy"]:
            continue
        g = r["Target_Gene"].strip().upper()
        if not g:
            continue
        a = agg[g]
        a["n"] += 1
        if r["InChI_Key"].strip():
            a["chemo"].add(r["InChI_Key"].strip())          # distinct compound = chemotype proxy
        if r["War_head"].strip():
            a["warheads"].add(r["War_head"].strip())
        site = r["Site"].strip()
        if site and site != "NT":
            a["residues"].add(site)
        if "kinact/ki" in r["Activity_type"].strip().lower():
            v = _num(r["Value"])
            if v is not None and (a["kk"] is None or v > a["kk"]):
                a["kk"] = v                                  # best (most efficient) covalent kinetics on record
    out = {g: {"n_inhibitors": a["n"], "n_chemotypes": len(a["chemo"]),
               "warheads": sorted(a["warheads"]), "residues": sorted(a["residues"]),
               "cys_sites": sorted(s for s in a["residues"] if s.upper().startswith("CYS")),
               "best_kinact_ki": a["kk"]}
           for g, a in agg.items()}
    json.dump(out, open(CACHE, "w"))
    return out


def lookup(gene):
    """Covalent prior-art for one gene, or None if CovalentInDB has no known covalent chemistry for it
    (None is NOT proof of openness — the cutoff is Jul-2024; pair with the clinical leg + competition check)."""
    return index().get((gene or "").upper())


if __name__ == "__main__":
    import sys
    if "--selftest" in sys.argv:
        if not CSV.exists():
            print("covindb selftest SKIPPED: CovInDB_All.csv not present (gitignored raw data)")
            sys.exit()
        idx = index()
        egfr = idx.get("EGFR")
        assert egfr and any("797" in s for s in egfr["cys_sites"]), f"EGFR CYS-797 expected: {egfr}"
        assert "WRN" not in idx, "WRN must be absent (verified recency gap — the TRAP case)"
        print(f"covindb self-check OK: {len(idx)} human targets; EGFR {egfr['n_inhibitors']} inh, sites {egfr['cys_sites'][:3]}")
    else:
        idx = index()
        print(f"{len(idx)} human covalent-prior-art targets -> {CACHE.name}")
        for g in ["EGFR", "BTK", "KRAS", "MALT1", "PRMT5", "WRN", "KEAP1"]:
            c = idx.get(g)
            print(f"  {g:7} " + (f"{c['n_inhibitors']:4} inh · {c['n_chemotypes']} chemo · {c['warheads'][:2]} · "
                                 f"cys {c['cys_sites'][:2]} · kk={c['best_kinact_ki']}" if c else "— none"))
