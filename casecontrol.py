"""CASE-CONTROL validation of patent-emergence as a target predictor (panel-prescribed; sidesteps the
~0% base rate of a random-universe test by using KNOWN emergences as positives).

POSITIVES: genes that became NOVEL drug targets via first-in-class clinical entry/approval 2018-2024
  (compiled by clinical milestone, NOT patent/pub buzz -> no leakage into the signal). In positives.json.
CONTROLS: genes from the random 3000-universe that did NOT emerge, matched to positives on patent LEVEL.
SIGNAL: patent-emergence = (patents 2015-2017 + 1)/(patents 2012-2014 + 1)  [strictly pre-2018 -> no leakage].

================= PRE-REGISTRATION (frozen before computing AUC) =================
PRIMARY: AUC of patent-emergence separating positives from p_recent-MATCHED controls (k=5 nearest by 2015-17
  patent level), one-sided Mann-Whitney p. SURVIVES iff AUC > 0.5 AND p < 0.05.
SECONDARY (context, not gating): AUC unmatched; AUC of raw recent patent COUNT (does emergence beat level?);
  later the same for grant- and literature-emergence on the same positives/controls.
INTERPRETATION GUARDRAIL: if matched-AUC ~ unmatched-AUC ~ recent-count-AUC, the "signal" is just patent LEVEL,
  not emergence. Emergence earns its keep only if matched-AUC(emergence) > AUC(recent count).
=================================================================================
"""
import bisect
import csv
import json
import math
from pathlib import Path

HERE = Path(__file__).parent
BASE_YEARS = range(2012, 2015)
RECENT_YEARS = range(2015, 2018)


def load_counts(path):
    d = {}
    for row in csv.DictReader(open(path)):
        d.setdefault(row["gene"], {})[int(row["yr"])] = int(row["n"])
    return d


def emergence(y):
    base = sum(y.get(yr, 0) for yr in BASE_YEARS)
    recent = sum(y.get(yr, 0) for yr in RECENT_YEARS)
    return recent, (recent + 1) / (base + 1)


def auc(pos, neg):                                  # P(pos > neg), ties=0.5  ==  Mann-Whitney U / (n_pos*n_neg)
    negs = sorted(neg)
    n = len(negs)
    tot = 0.0
    for p in pos:
        lo = bisect.bisect_left(negs, p)
        hi = bisect.bisect_right(negs, p)
        tot += lo + 0.5 * (hi - lo)
    return tot / (len(pos) * n) if pos and n else 0.5


def mw_p(pos, neg):                                 # one-sided (pos > neg) Mann-Whitney, normal approx
    n1, n2 = len(pos), len(neg)
    U = auc(pos, neg) * n1 * n2
    mu = n1 * n2 / 2
    sd = math.sqrt(n1 * n2 * (n1 + n2 + 1) / 12) or 1
    z = (U - mu) / sd
    return 0.5 * math.erfc(z / math.sqrt(2))


def match_controls(positives, control_pool, k=5):   # k nearest controls per positive by p_recent level, no replacement
    used, chosen = set(), []
    pool = sorted(control_pool, key=lambda c: c["p_recent"])
    levels = [c["p_recent"] for c in pool]
    for p in positives:
        idx = bisect.bisect_left(levels, p["p_recent"])
        picks, lo, hi = [], idx - 1, idx
        while len(picks) < k and (lo >= 0 or hi < len(pool)):
            cand = None
            if hi >= len(pool):
                cand = lo
                lo -= 1
            elif lo < 0:
                cand = hi
                hi += 1
            elif abs(pool[lo]["p_recent"] - p["p_recent"]) <= abs(pool[hi]["p_recent"] - p["p_recent"]):
                cand = lo
                lo -= 1
            else:
                cand = hi
                hi += 1
            if pool[cand]["gene"] not in used:
                used.add(pool[cand]["gene"])
                picks.append(pool[cand])
        chosen += picks
    return chosen


def run():
    pc = load_counts(HERE / "patent_counts.csv")
    positives_raw = json.loads((HERE / "positives.json").read_text())
    pos_syms = {p["symbol"] for p in positives_raw}
    universe = set(json.loads((HERE / "universe.json").read_text()))

    pos = []
    for s in pos_syms:
        if s in pc:
            r, e = emergence(pc[s])
            pos.append({"gene": s, "p_recent": r, "emergence": e})
    controls = []
    for g in universe:
        if g in pos_syms:
            continue
        r, e = emergence(pc.get(g, {}))
        controls.append({"gene": g, "p_recent": r, "emergence": e})

    matched = match_controls(pos, controls, k=5)
    pe = [p["emergence"] for p in pos]
    ce_all = [c["emergence"] for c in controls]
    ce_m = [c["emergence"] for c in matched]
    pr = [p["p_recent"] for p in pos]
    cr_all = [c["p_recent"] for c in controls]

    print(f"positives with patent data: {len(pos)}/{len(pos_syms)}  | controls pool: {len(controls)}  | matched controls: {len(matched)}")
    print(f"median patent-emergence  positives={_med(pe):.2f}  matched-controls={_med(ce_m):.2f}  all-controls={_med(ce_all):.2f}")
    a_un = auc(pe, ce_all); a_m = auc(pe, ce_m); a_lvl = auc(pr, cr_all)
    print(f"AUC patent-EMERGENCE vs matched controls = {a_m:.3f}  (one-sided MW p = {mw_p(pe, ce_m):.4f})   [PRIMARY]")
    print(f"AUC patent-EMERGENCE vs ALL controls      = {a_un:.3f}  (p = {mw_p(pe, ce_all):.4f})")
    print(f"AUC raw recent patent COUNT vs all        = {a_lvl:.3f}   <- emergence must beat THIS to matter")
    verdict = "SURVIVES" if (a_m > 0.5 and mw_p(pe, ce_m) < 0.05) else "FAILS"
    beats_level = "and beats level-AUC" if a_m > a_lvl else "but does NOT beat level-AUC (signal is just patent count)"
    print(f"PRE-REGISTERED VERDICT: patent-emergence {verdict} ({beats_level})")
    print("\npositives ranked by patent-emergence (signal as-of-2017 -> these DID emerge 2018-24):")
    for p in sorted(pos, key=lambda x: x["emergence"], reverse=True)[:20]:
        print(f"  {p['gene']:9} emergence={p['emergence']:6.1f}  recent(2015-17)={p['p_recent']:4}")


def _med(xs):
    import statistics
    return statistics.median(xs) if xs else 0.0


if __name__ == "__main__":
    run()
