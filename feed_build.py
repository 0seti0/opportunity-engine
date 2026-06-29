"""Leg 1 — the Day-1 covalent-patent screen. Ingest by TARGET + pharma IPC (not by warhead keyword),
pull WO/EP full-text, and decide covalency from the warhead NOMENCLATURE in the text — robustly, without
depending on any "Example N" scaffold (real patents label syntheses "Step 7:", "Compound 12", etc.).

Signals (learned against the live OPS API — see DESIGN.md):
  - TRIGGER: count of syntax-bound warhead names in the description (acrylamide/prop-2-enamide/but-2-ynamide/
    chloroacetamide/…), reusing patent_radar's nomenclature regex (its exclusion list kills 'polyacrylamide gel'
    polymer prose). Clean separation observed: a covalent series scores ~30+, a non-covalent patent scores 0.
  - EVIDENCE: pull the warhead-bearing compound NAMES, repair OPS text artifacts ('/V' -> italic 'N'),
    resolve with local OPSIN, confirm with warhead SMARTS -> the actual lead structures.
  - HANDOFF: when the text says covalent (mentions high) but no name resolves, the warheads are DRAWN, not
    written -> needs_ocsr=True routes the patent to the image OCSR path (ocsr.py / option 3).

Why this shape: keyword-searching the warhead returns polymer art; OPS sort is CN-heavy with no CN full-text
(post-filter to WO/EP: 25% -> 100% coverage); OPS description is text-only IUPAC names (no SMILES/images).
"""
import html
import re
import time

import ops
import patent_radar as pr

IPC = "A61P35"                 # antineoplastic; widen per target class (A61P29/37 for immuno/inflammation)
MIN_MENTIONS = 3               # >= this many syntax-bound warhead names in the text => covalent series

# covalent-warhead NAME suffixes (the forms that appear inside a compound's IUPAC name, not reagents/prose)
_WARHEAD_NAME = re.compile(r"prop-2-enamide|acrylamide|but-2-ynamide|propiolamide|cyanoacrylamide|"
                           r"prop-2-enenitrile|maleimid", re.I)


def _clean_name(p):
    """Repair OPS full-text artifacts and strip a leading synthesis label / trailing analytical data."""
    p = p.replace("/V", "N").replace("/v", "N")              # OPS renders italic 'N' (locant prefix) as '/V'
    p = re.sub(r"^(?:step\s+\d+[:.]?|example\s+\d+[A-Za-z]*[:.]?|compound\s+\d+[:.]?|title compound[,:]?)\s*",
               "", p, flags=re.I)
    return re.split(r"\s+(?:1H ?NMR|13C|LCMS|LC-MS|MS ?\(|ESI|δ|m/?z|HRMS|Yield|as a (?:white|pale|yellow|off|colou?r))",
                    p)[0].strip()


def warhead_leads(desc_xml):
    """Cleaned compound names carrying a warhead suffix — the covalent leads named in the description text."""
    paras = [re.sub(r"\s+", " ", html.unescape(re.sub(r"<(\d+)>", r"\1", p))).strip()  # <3> subscript -> 3
             for p in re.findall(r"<p>(.*?)</p>", desc_xml, re.S)]
    leads = [_clean_name(p) for p in paras
             if 12 < len(p) < 320 and not p.startswith("[") and _WARHEAD_NAME.search(p)]
    return list(dict.fromkeys(leads))


def screen(target, ipc=IPC, max_pubs=25, pause=1.0):
    """Screen recent WO/EP inhibitor patents for <target>; return covalent-patent events."""
    pubs = [p for p in ops.search(f"ic={ipc} and ti=inhibitor and txt={target}", n=max_pubs)
            if p[:2] in ("WO", "EP")]                         # CN/US members 404 on OPS full-text
    events = []
    for pn in pubs:
        time.sleep(pause)                                    # gentle on the OPS retrieval bucket (50/min)
        desc = ops.fulltext(pn, "description")
        if not desc:
            continue
        desc_text = re.sub(r"<[^>]+>", " ", desc)
        if target.upper() not in desc_text.upper():
            continue                                         # cross-target bleed guard (txt= matches anywhere)
        mentions = pr.warhead_in_claims(desc_text)           # syntax-bound warhead nomenclature, polymer-excluded
        leads = warhead_leads(desc)
        smis = pr.resolve_batch(leads)
        confirmed = [{"name": n[:90], "smiles": s, "warhead": pr.has_warhead(s)}
                     for n, s in zip(leads, smis) if s and pr.has_warhead(s)]
        claims = ops.fulltext(pn, "claims")
        claim_wh = pr.warhead_in_claims(claims) if claims else []
        if len(mentions) >= MIN_MENTIONS or confirmed or claim_wh:
            events.append({"pn": pn, "target": target.upper(),
                           "warhead_mentions": len(mentions), "n_warhead_names": len(leads),
                           "n_confirmed": len(confirmed), "named_in_claims": claim_wh[:2],
                           # text says covalent but nothing resolved -> warheads are DRAWN -> OCSR (ocsr.py)
                           "needs_ocsr": len(mentions) >= MIN_MENTIONS and not confirmed,
                           "leads": confirmed[:5]})
    return events


def verify(pn, target):
    """AUTHORITATIVE per-patent attribution check. feed_build is the authority (target-centric + composition-
    of-matter + warhead); attribute.py's text-mined grade is only a CANDIDATE that this confirms. Fetches the
    OPS full-text and applies the screen's own covalency + on-target test. `pn` may be a SureChEMBL id like
    'EP-4067347-B1'; only WO/EP carry OPS full-text (others return verified=False, reason given)."""
    epodoc = re.sub(r"^([A-Z]{2}\d+).*", r"\1", (pn or "").replace("-", ""))   # 'EP-4067347-B1' -> 'EP4067347'
    desc = ops.fulltext(epodoc, "description")
    if not desc:
        return {"pn": pn, "target": target.upper(), "verified": False, "reason": "no OPS full-text (WO/EP only)"}
    text = re.sub(r"<[^>]+>", " ", desc)
    covalent = len(pr.warhead_in_claims(text)) >= MIN_MENTIONS
    on_target = target.upper() in text.upper()
    return {"pn": pn, "target": target.upper(), "verified": bool(covalent and on_target),
            "covalent": covalent, "target_in_text": on_target}


def ocsr_rescue(pn, frac_range=(0.4, 0.95), max_pages=8):
    """For a `needs_ocsr` patent (covalent by text, but no compound NAME resolved -> warheads are DRAWN):
    recover them from the description page images. Returns the warhead-bearing structures DECIMER read.
    Heavy (DECIMER/TF) and OPS-image-call-heavy, so it's a deliberate second pass, not part of screen()."""
    import ocsr
    inv = ops.images_inventory(pn)
    if not inv or "DESCRIPTION" not in inv["sections"]:
        return []
    window = list(range(inv["sections"]["DESCRIPTION"], inv["sections"].get("CLAIMS", inv["pages"])))
    lo, hi = (int(len(window) * f) for f in frac_range)        # examples sit in the back of the description
    pages = window[lo:hi][:: max(1, (hi - lo) // max_pages)][:max_pages]
    out = []
    for p in pages:
        time.sleep(0.5)
        pdf = ops.page_pdf(inv["link"], p)
        if not pdf:
            continue
        for s in ocsr.page_structures(ocsr.pdf_to_png(pdf)):
            if pr.has_warhead(s):
                out.append({"page": p, "smiles": s, "warhead": pr.has_warhead(s)})
    return out


if __name__ == "__main__":
    import json
    import sys

    # offline parser self-check (deterministic): /V->N repair, label strip, trailing-data trim, prose exclusion
    fixture = ("<p>Step 7: (E)-/V-phenyl-3-(pyrrolidin-2-yl)prop-2-enamide as a white solid (0.03 g)</p>"
               "<p>[0050] To a stirred solution of the title compound in THF...</p>")
    leads = warhead_leads(fixture)
    assert leads == ["(E)-N-phenyl-3-(pyrrolidin-2-yl)prop-2-enamide"], leads
    assert re.sub(r"^([A-Z]{2}\d+).*", r"\1", "EP-4067347-B1".replace("-", "")) == "EP4067347"  # verify() pn norm
    print("parser self-check OK:", leads)

    target = sys.argv[1] if len(sys.argv) > 1 else "EGFR"
    print(f"\nscreening {target} (live OPS) ...")
    evs = screen(target)
    for e in evs:
        flag = " [DRAWN→OCSR]" if e["needs_ocsr"] else ""
        print(f"  ⚡ {e['pn']} {e['target']}  mentions={e['warhead_mentions']}  "
              f"confirmed-leads={e['n_confirmed']}/{e['n_warhead_names']}  claims={e['named_in_claims']}{flag}")
        for ld in e["leads"][:2]:
            print(f"       {ld['warhead']} {ld['smiles']}")
    print(f"{len(evs)} covalent event(s)")
