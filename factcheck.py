"""Deterministic structural-fact verification layer — NO LLM, so it cannot hallucinate (the PRKCQ
'C303 is not a cysteine' error was an LLM agent asserting a sequence fact; here it is checked in code).

For each opportunity's claimed covalent residue: (1) confirm it IS that amino acid in the UniProt
sequence; (2) align to paralogs and read the aligned residue -> selectivity = unique-handle (paralogs
differ) vs conserved-liability (paralogs share the Cys) vs partial. Pure UniProt + Biopython alignment.
"""
import json
import urllib.request
from pathlib import Path

from Bio import Align
from Bio.Align import substitution_matrices

HERE = Path(__file__).parent
CACHE = HERE / ".uniprot_cache.json"
_seqs = json.loads(CACHE.read_text()) if CACHE.exists() else {}

# target -> (uniprot, [claimed Cys positions], {paralog: uniprot})
TARGETS = {
    "SLC1A5": ("Q15758", [467], {"SLC1A4": "P43007", "SLC1A1": "P43005", "SLC1A2": "P43004",
                                 "SLC1A3": "P43003", "SLC1A6": "P48664", "SLC1A7": "O00341"}),
    "PPARD":  ("Q03181", [249], {"PPARA": "Q07869", "PPARG": "P37231"}),
    "USP35":  ("Q9P2H5", [450], {"USP30": "Q70CQ3", "USP7": "Q93009"}),
    "IRF4":   ("Q15306", [99, 194, 250], {"IRF8": "Q02556", "IRF1": "P10914", "IRF3": "Q14653"}),
    "YES1":   ("P07947", [287], {"SRC": "P12931", "FGR": "P09769", "FYN": "P06241", "LCK": "P06239"}),
    "RAF1":   ("P04049", [424], {"BRAF": "P15056", "ARAF": "P10398"}),
    "AKR1B1": ("P15121", [299], {"AKR1B10": "O60218"}),   # lit 'Cys298' is mature numbering = UniProt 299
    "PRKCG":  ("P05129", [516], {"PRKCA": "P17252", "PRKCB": "P05771"}),
    "PRKCQ":  ("Q04759", [303, 322, 424], {"PRKCD": "Q05655", "PRKCE": "Q02156", "PRKCH": "P24723"}),
}

_aligner = Align.PairwiseAligner(mode="global", open_gap_score=-11, extend_gap_score=-1)
_aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")


def seq(acc):
    if acc not in _seqs:
        url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
        txt = urllib.request.urlopen(url, timeout=45).read().decode()
        _seqs[acc] = "".join(l.strip() for l in txt.splitlines() if not l.startswith(">"))
        CACHE.write_text(json.dumps(_seqs))
    return _seqs[acc]


def aligned_residue(tseq, pos, pseq):           # residue in pseq aligned to tseq position `pos` (1-based)
    aln = _aligner.align(tseq, pseq)[0]
    for (ts, te), (ps, pe) in zip(*aln.aligned):
        if ts <= pos - 1 < te:                  # pos-1 = 0-based index inside this aligned block
            return pseq[ps + (pos - 1 - ts)]
    return "-"                                  # falls in a gap


def run():
    out = {}
    for sym, (acc, positions, paras) in TARGETS.items():
        ts = seq(acc)
        rec = {"uniprot": acc, "residues": [], "selectivity": None}
        sel_votes = []
        for p in positions:
            is_cys = p <= len(ts) and ts[p - 1] == "C"
            cons = {}
            for pn, pacc in paras.items():
                try:
                    cons[pn] = aligned_residue(ts, p, seq(pacc))
                except Exception as e:
                    cons[pn] = f"err"
            ncys = sum(1 for r in cons.values() if r == "C")
            verdict = ("unique-handle" if ncys == 0 else
                       "conserved-liability" if ncys == len(cons) and cons else "partial")
            sel_votes.append(verdict)
            rec["residues"].append({"pos": p, "is_cysteine": is_cys, "paralog_residues": cons, "selectivity": verdict})
        # best (most selective) handle drives the target-level selectivity call
        rec["selectivity"] = min(sel_votes, key=lambda v: {"unique-handle": 0, "partial": 1, "conserved-liability": 2}[v]) if sel_votes else "n/a"
        out[sym] = rec

    json.dump(out, open(HERE / "factcheck.json", "w"), indent=1)
    print(f"{'target':8}{'residue checks':<40}{'selectivity':>20}")
    for sym, rec in out.items():
        checks = "  ".join(f"C{r['pos']}={'✓' if r['is_cysteine'] else '✗NOT-CYS'}"
                           + "(" + "".join(rec_v for rec_v in [{k: v for k, v in r['paralog_residues'].items()}][0].values()) + ")"
                           for r in rec["residues"]) or "(no positions)"
        print(f"{sym:8}{checks:<40}{rec['selectivity']:>20}")
    print("\nsaved factcheck.json")
    # the headline corrections
    print("\n-- key checks --")
    pq = out["PRKCQ"]["residues"]
    print(f"PRKCQ: C303={pq[0]['is_cysteine']} C322={pq[1]['is_cysteine']} C424={pq[2]['is_cysteine']}  (agent claimed these are NOT cysteines)")
    s = out["SLC1A5"]["residues"][0]
    print(f"SLC1A5 C467: is_cys={s['is_cysteine']} paralogs={s['paralog_residues']} -> {s['selectivity']}  (the #1 claim)")


if __name__ == "__main__":
    run()
