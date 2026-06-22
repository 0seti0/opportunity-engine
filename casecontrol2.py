"""Case-control, ALL signals — does ANY public signal (patent/grant/literature, emergence OR level) separate
genes that became novel drug targets from patent-LEVEL-matched controls? (the fair 'among candidates' test).

Matching controls on patent level neutralizes the trivial 'targets have patents' effect, so a signal that
still separates positives there is a genuine edge BEYOND raw patent volume. All signal windows are pre-2018
(base 2012-2014, recent 2015-2017) -> no leakage. Network: grant (NIH RePORTER) + literature (Europe PMC)
year counts for positives + matched controls only (cached).
"""
import json
import time
from pathlib import Path

import backtest as bt
import casecontrol as cc

HERE = Path(__file__).parent
BASE = (2012, 2013, 2014)
RECENT = (2015, 2016, 2017)


def _safe(fn, *a):                                  # retry hard with backoff; None if it never succeeds
    for attempt in range(6):
        try:
            return fn(*a)
        except Exception:
            if attempt < 5:
                time.sleep(4 * (attempt + 1)); continue
            return None


def add_grant_lit(rows):
    for i, row in enumerate(rows):
        g = row["gene"]
        gy = [_safe(bt.reporter_year, g, y) for y in BASE + RECENT]
        ly = [_safe(bt.epmc_year, g, y) for y in BASE + RECENT]
        if None not in gy:
            gb, gr = sum(gy[:3]), sum(gy[3:])
            row["g_recent"], row["g_emerg"] = gr, (gr + 1) / (gb + 1)
        if None not in ly:
            lb, lr = sum(ly[:3]), sum(ly[3:])
            row["l_recent"], row["l_emerg"] = lr, (lr + 1) / (lb + 1)
        time.sleep(0.15)                            # gentle throttle so EPMC/RePORTER don't rate-limit the burst
        if i % 10 == 9:
            bt.CACHE.write_text(json.dumps(bt._cache))
    bt.CACHE.write_text(json.dumps(bt._cache))


def run():
    pc = cc.load_counts(HERE / "patent_counts.csv")
    pos_syms = {p["symbol"] for p in json.loads((HERE / "positives.json").read_text())}
    universe = set(json.loads((HERE / "universe.json").read_text()))

    def rec(g):
        r, e = cc.emergence(pc.get(g, {}))
        return {"gene": g, "p_recent": r, "p_emerg": e}

    pos = [rec(g) for g in pos_syms if g in pc]
    controls = [rec(g) for g in universe if g not in pos_syms]
    matched = cc.match_controls(pos, controls, k=5)

    add_grant_lit(pos)
    add_grant_lit(matched)

    def line(name, key):
        pv = [p[key] for p in pos if key in p]
        cv = [c[key] for c in matched if key in c]
        if not pv or not cv:
            print(f"  {name:24} (insufficient data)"); return
        a = cc.auc(pv, cv); pval = cc.mw_p(pv, cv)
        cav = [c[key] for c in controls if key in c]
        extra = f"   AUC_vs_all={cc.auc(pv, cav):.3f}" if len(cav) > 100 else ""
        flag = "  <== beats chance" if (a > 0.5 and pval < 0.05) else ""
        print(f"  {name:24} AUC_matched={a:.3f} (p={pval:.3f}, n_pos={len(pv)}){extra}{flag}")

    print(f"positives={len(pos)}  matched controls={len(matched)}  (matched on patent level)\n")
    print("Signal separation, positives vs patent-LEVEL-matched controls (the fair test):")
    for nm, k in [("patent emergence", "p_emerg"), ("patent level", "p_recent"),
                  ("grant emergence", "g_emerg"), ("grant level", "g_recent"),
                  ("literature emergence", "l_emerg"), ("literature level", "l_recent")]:
        line(nm, k)

    # combined level: rank-sum of patent + grant level over rows that have both
    keys = ["p_recent", "g_recent"]
    both = [r for r in pos + matched if all(k in r for k in keys)]
    for k in keys:
        order = sorted(range(len(both)), key=lambda i: both[i][k])
        for rank_i, i in enumerate(order):
            both[i][k + "_r"] = rank_i
    posset = {p["gene"] for p in pos}
    pcomb = [r["p_recent_r"] + r["g_recent_r"] for r in both if r["gene"] in posset]
    ccomb = [r["p_recent_r"] + r["g_recent_r"] for r in both if r["gene"] not in posset]
    if pcomb and ccomb:
        print(f"\n  combined patent+grant LEVEL (rank-sum)  AUC_matched={cc.auc(pcomb, ccomb):.3f}")
    print("\nINTERPRETATION: any AUC_matched >> 0.5 = a signal that beats raw patent volume among candidates.")


if __name__ == "__main__":
    run()
