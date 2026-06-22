"""Drug resolver -> DrugResolution (resolved / ambiguous / unresolved), mirroring the gene resolver.
Dual source (GtoPdb + ChEMBL); collapses salt/variant names to one drug via the ChEMBL parent
(metformin == metformin hydrochloride); target symbols canonicalize to HGNC ids through Phase 1."""
import json
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
from typing import Optional

from identity import Outcome, build_index, load_rows, resolve  # reuse Phase 1

GTP = "https://www.guidetopharmacology.org/services/"
CHEMBL = "https://www.ebi.ac.uk/chembl/api/data/"


@dataclass
class DrugResolution:
    outcome: Outcome
    drug_id: Optional[str] = None                    # canonical drug (ChEMBL parent) when RESOLVED
    targets: list = field(default_factory=list)      # target HGNC ids when RESOLVED
    candidates: list = field(default_factory=list)   # candidate drug ids when AMBIGUOUS


def _get(url):
    return json.load(urllib.request.urlopen(
        urllib.request.Request(url, headers={"User-Agent": "x", "Accept": "application/json"}), timeout=30))


def _get_or_empty(url, empty):   # 404 = deterministic "no match" -> empty; other errors re-raise
    try:
        return _get(url)
    except urllib.error.HTTPError as e:
        if e.code == 404:
            return empty
        raise


def _gtp(path, **p):
    return _get_or_empty(GTP + path + ("?" + urllib.parse.urlencode(p) if p else ""), [])


def _chembl(path, **p):
    p["format"] = "json"
    return _get_or_empty(CHEMBL + path + "?" + urllib.parse.urlencode(p), {})


def _parent(m):                  # a salt's parent drug id, or the molecule itself
    return (m.get("molecule_hierarchy") or {}).get("parent_chembl_id") or m.get("molecule_chembl_id")


def _chembl_resolve(name):
    """-> (distinct parent drugs the exact name maps to, single-protein target gene symbols for that drug).
    Gathers mechanisms by PARENT, so name variants (metformin vs metformin hydrochloride) give the SAME result.
    ponytail: no throttle/cache yet -- add when batching thousands of drugs (Phase 0 has the pattern)."""
    mols = _chembl("molecule/search", q=name).get("molecules", [])
    if not mols:
        return set(), set()
    exact = [m for m in mols if (m.get("pref_name") or "").upper() == name.upper()] or [mols[0]]
    parents = {_parent(m) for m in exact}
    if len(parents) != 1:
        return parents, set()                              # >1 distinct drug -> caller flags AMBIGUOUS
    genes = set()
    for mech in _chembl("mechanism", parent_molecule_chembl_id=next(iter(parents))).get("mechanisms", []):
        tid = mech.get("target_chembl_id")
        if not tid:
            continue
        t = _chembl(f"target/{tid}")
        if t.get("target_type") != "SINGLE PROTEIN":       # collapse: a complex/family/fusion isn't one gene-level target
            continue
        for comp in t.get("target_components", []):
            for syn in comp.get("target_component_synonyms", []):
                if syn.get("syn_type") == "GENE_SYMBOL":
                    genes.add(syn["component_synonym"])
    return parents, genes


def _gtp_targets(name):
    ligs = _gtp("ligands", name=name)
    if not ligs:
        return set()
    genes = set()
    for it in _gtp(f"ligands/{ligs[0]['ligandId']}/interactions"):
        for row in _gtp(f"targets/{it['targetId']}/geneProteinInformation", species="Human"):
            if row.get("geneSymbol"):
                genes.add(row["geneSymbol"])
    return genes


SALTS = ("hydrochloride", "dihydrochloride", "hydrobromide", "sulfate", "sulphate", "sodium",
         "calcium", "potassium", "mesylate", "maleate", "tartrate", "acetate", "citrate",
         "succinate", "fumarate", "phosphate", "besylate", "monohydrate", "dihydrate", "anhydrous")


def _base_name(name):   # ponytail: strip trailing salt/form words so "metformin hydrochloride" == "metformin"; may mis-strip rare names
    parts = name.strip().lower().split()
    while len(parts) > 1 and parts[-1] in SALTS:
        parts.pop()
    return " ".join(parts)


def resolve_drug(name, gene_index):
    name = _base_name(name)                          # collapse salt/form variants across BOTH sources
    parents, c_genes = _chembl_resolve(name)
    if len(parents) > 1:                             # name maps to several distinct drugs
        return DrugResolution(Outcome.AMBIGUOUS, candidates=sorted(parents))
    symbols = c_genes or _gtp_targets(name)          # GtoPdb only as fallback when ChEMBL has no target (faster; keeps code-name coverage)
    if not symbols:                              # resolved means we found a target; no target (incl. fuzzy junk) -> unresolved
        return DrugResolution(Outcome.UNRESOLVED)
    targets = set()
    for s in symbols:                                # canonicalize target symbols -> HGNC ids via Phase 1
        hid = resolve(s, gene_index).hgnc_id
        if hid:
            targets.add(hid)
    drug_id = sorted(parents)[0] if parents else None
    return DrugResolution(Outcome.RESOLVED, drug_id=drug_id, targets=sorted(targets))


if __name__ == "__main__":
    idx = build_index(load_rows())
    _cache = {}
    def R(n):
        if n not in _cache:
            _cache[n] = resolve_drug(n, idx)
        return _cache[n]
    assert R("trastuzumab").outcome is Outcome.RESOLVED and "HGNC:3430" in R("trastuzumab").targets  # single protein kept
    assert R("totally-not-a-drug").outcome is Outcome.UNRESOLVED
    assert len(R("metformin").targets) < 5                                  # Complex-I family collapsed (was 51)
    assert R("imatinib").targets and R("imatinib").targets == R("imatinib mesylate").targets  # salt collapse, single proteins kept
    for n in ["trastuzumab", "imatinib", "metformin", "totally-not-a-drug"]:
        r = R(n)
        print(f"  {n:22} -> {r.outcome.value:10} {r.targets}")
    print("drug resolver checks green")
