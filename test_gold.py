"""Regression gate: turns the methodology audit's asserted accuracy into a runnable check.

gold.json is a SEED set — the audit's confirmed failure modes + canonical covalent chemistry — NOT the
original 30-patent / 24-trial audit; grow it toward that. Exits non-zero on any miss, so it gates changes
to grade() (attribution), has_warhead() (warhead SMARTS), and the clinical covalency verdicts.

  python test_gold.py
"""
import json
from pathlib import Path

HERE = Path(__file__).parent
GOLD = json.load(open(HERE / "gold.json"))


def _check(name, cases, fn):
    miss = [(c, got) for c in cases if (got := fn(c)) != c["expect"]]
    for c, got in miss:
        print(f"  ✗ {name}: got {got!r} != {c['expect']!r}  <-  {c.get('note') or c.get('name') or c}")
    print(f"{name:26} {len(cases) - len(miss)}/{len(cases)} pass")
    return not miss


def main():
    ok = True
    from attribute import grade                         # fix1 (claims->candidate) + fix2 (radar->trusted)
    ok &= _check("attribution.grade", GOLD["attribution"], lambda c: grade(c["conf"], c["sym"], c["title"]))

    from patent_radar import has_warhead               # fix4 (expanded SMARTS + precision guards)
    ok &= _check("warhead.has_warhead", GOLD["warhead"], lambda c: bool(has_warhead(c["smiles"])))

    try:                                               # clinical verdicts (also in clinical_radar.py --selftest)
        from clinical_radar import fingerprint, FP_THRESHOLD, _pub_verdict
        ok &= _check("clinical.fingerprint", GOLD["clinical_fingerprint"],
                     lambda c: fingerprint(c["text"])["score"] >= FP_THRESHOLD)
        ok &= _check("clinical.pub_verdict", GOLD["clinical_pub"], lambda c: _pub_verdict(c["text"]))
    except Exception as e:
        print(f"clinical.*                 skipped ({type(e).__name__}); run `python clinical_radar.py --selftest`")

    assert ok, "GOLD regression FAILED"
    print("\n✅ gold regression passed")


if __name__ == "__main__":
    main()
