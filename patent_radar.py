"""Step 0+1 of the covalent early-warning radar — the Day-1 patent screen (no SureChEMBL wait).

  independent_claims(xml) -> the independent claims (where the CORE invention lives; we search only here)
  warhead_in_claims(text) -> syntax-bound warhead regex + exclusion list  (catches NAMED warheads)
  resolve_to_smiles(rep)  -> a compound (SMILES/InChI | IUPAC name | image) -> SMILES, preferring the most
                             reliable representation (RDKit | OPSIN | DECIMER)
  has_warhead(smiles)     -> covalent-warhead SMARTS                       (catches DRAWN warheads)

Runs on raw patent XML at publication. SureChEMBL becomes a later cross-check, not the source.
"""
import re
import xml.etree.ElementTree as ET

from rdkit import Chem

# --- covalent-warhead SMARTS (acyclic constraints exclude steroid-enone / fumarate false positives) -------
WARHEAD_SMARTS = {
    # Michael-acceptor / alkylating amides (acrylamide subsumes cyano- & fluoro-acrylamide — same C=C-C(=O)N core)
    "acrylamide":        "[CX3;!R]=[CX3;!R]C(=O)[NX3]",
    "haloacetamide":     "[Cl,Br,I]C[CX3](=O)[NX3]",        # chloro/bromo/iodo-acetamide (was chloro-only)
    "propiolamide":      "[CX2]#[CX2]C(=O)[NX3]",
    "vinylsulfonamide":  "[CX3;!R]=[CX3;!R]S(=O)(=O)[NX3]",
    "maleimide":         "O=C1C=CC(=O)N1",                  # thiol-Michael (also in the text-name regex -> paths agree)
    # clinically-validated covalent classes the original 4-pattern set left structurally invisible (audit):
    "sulfonyl_fluoride": "[#16X4](=O)(=O)[F]",              # SuFEx warhead
    "epoxide":           "[CX4]1[OX2][CX4]1",               # oxirane / epoxyketone (carfilzomib, proteasome)
    "aziridine":         "[CX4]1[NX3][CX4]1",
    "beta_lactam":       "[NX3]1[CX3](=O)[CX4][CX4]1",
    "boronic_acid":      "[#5]([OX2])[OX2]",                # reversible-covalent (bortezomib/ixazomib, proteasome)
}
# Deliberately NOT added — precision killers in the SureChEMBL-wide structure scan, where a greedy pattern over
# ALL exemplified compounds (>=MIN_CPDS) would mass-false-positive: generic nitrile [CX2]#[NX1] (most drug
# nitriles are non-covalent), bare aldehyde (ubiquitous/unstable), activated-heteroaryl SNAr halide (an aryl-Cl
# synthetic handle, not a warhead). Add context-gated if a real miss surfaces.  # ponytail: precision ceiling
_SMARTS = {k: Chem.MolFromSmarts(v) for k, v in WARHEAD_SMARTS.items()}

# --- claims-text tripwire: warhead NAMES (common + IUPAC), only when they read as chemical nomenclature ----
_NAME = re.compile(r"(?:acryl|prop-?2-?en|2-?propen)amide|(?:propiol|but-?2-?yn)amide|"
                   r"chloroacetamide|vinyl[\s-]?sulfon(?:amide|e)|epoxide|maleimid", re.I)
_EXCLUDE = re.compile(r"\b(gel|poly|polymer|excipient|column|purification|resin|hydrogel|coating)\b", re.I)


def warhead_in_claims(text):
    """Warhead names that read as nomenclature (preceded by a locant/bracket/hyphen), not prose."""
    hits = []
    for m in _NAME.finditer(text or ""):
        pre = text[max(0, m.start() - 2):m.start()]
        if re.search(r"[\]\)\-0-9]\s?$", pre):                 # chemical-syntax boundary just before the name
            ctx = text[max(0, m.start() - 28):m.end() + 12].strip()
            if not _EXCLUDE.search(ctx):
                hits.append(ctx)
    return hits


def independent_claims(xml_text):
    """Independent claims from patent XML. Handles two real shapes:
       ST.36     — one <claim> element per claim (dependent = <claim-ref> or 'claim N').
       EPO OPS   — one <claims>/<claim> blob whose <claim-text> lines are numbered '1.', '2.'.
    A claim is DEPENDENT if it back-references another claim."""
    root = ET.fromstring(xml_text)
    ln = lambda el: el.tag.split('}')[-1]
    claims = [el for el in root.iter() if ln(el) == "claim"]

    def independent(txt):                                   # 'of claim 1', 'any of claims 1-5' -> dependent
        return not re.search(r"\bclaims?\s+\d+", txt, re.I)

    out = []
    if len(claims) == 1 and sum(1 for c in claims[0].iter() if ln(c) == "claim-text") > 1:
        # OPS full-text: split the single blob by leading 'N.' numbering across its claim-text lines
        units, cur = [], []
        for line in ("".join(ct.itertext()).strip() for ct in claims[0].iter() if ln(ct) == "claim-text"):
            if re.match(r"^\d+\.\s", line) and cur:
                units.append(" ".join(cur)); cur = []
            cur.append(line)
        if cur:
            units.append(" ".join(cur))
        out = [re.sub(r"\s+", " ", u).strip() for u in units if independent(u)]
    else:
        for cl in claims:                                   # ST.36: each <claim> is one claim
            txt = re.sub(r"\s+", " ", "".join(cl.itertext())).strip()
            if not any(ln(c) == "claim-ref" for c in cl.iter()) and independent(txt):
                out.append(txt)
    return out


def _opsin(name):
    """IUPAC name -> SMILES via local OPSIN (py2opsin's bundled jar + a JRE). Deterministic, offline —
    no network, so a run is reproducible regardless of opsin.ch / cactus uptime or version drift."""
    from py2opsin import py2opsin
    return py2opsin(name) or None                              # '' (unparseable name) -> None


def resolve_to_smiles(rep, kind):
    """rep -> canonical SMILES. kind: smiles | inchi | iupac | image. Prefer the most reliable available."""
    if kind == "smiles":
        m = Chem.MolFromSmiles(rep)
    elif kind == "inchi":
        m = Chem.MolFromInchi(rep)
    elif kind == "iupac":
        s = _opsin(rep)
        m = Chem.MolFromSmiles(s) if s else None
    else:                                                       # image (PNG path) -> DECIMER OCSR (see ocsr.py)
        import ocsr
        return ocsr.image_to_smiles(rep)
    return Chem.MolToSmiles(m) if m else None


def resolve_batch(names):
    """IUPAC names -> canonical SMILES (None per failure), in ONE local-OPSIN/java call for the whole list."""
    from py2opsin import py2opsin
    if not names:
        return []
    raw = py2opsin(list(names))                               # list in -> list out, '' for an unparseable name
    out = []
    for s in (raw if isinstance(raw, list) else [raw]):
        m = Chem.MolFromSmiles(s) if s else None
        out.append(Chem.MolToSmiles(m) if m else None)
    return out


def has_warhead(smiles):
    m = Chem.MolFromSmiles(smiles) if smiles else None
    return [k for k, p in _SMARTS.items() if m and m.HasSubstructMatch(p)]


if __name__ == "__main__":
    print("=== 1. claims-text tripwire (nomenclature vs prose) ===")
    for t in ["wherein said compound is N-(3-chlorophenyl)prop-2-enamide and a pharmaceutically",
              "the core scaffold is substituted with a 2-chloroacetamide moiety at the R3 position",
              "the product was separated on a polyacrylamide gel for purification and analysis"]:
        print(f"   {str(warhead_in_claims(t)) or '[no warhead]':38} <- {t[:52]}")

    print("\n=== 2. independent vs dependent claims (search only the independent) ===")
    xml = ('<patent><claims>'
           '<claim id="c1"><claim-text>A compound of formula (I) bearing an acrylamide warhead, or a salt thereof.</claim-text></claim>'
           '<claim id="c2"><claim-text>The compound of <claim-ref idref="c1"/>, wherein R1 is methyl.</claim-text></claim>'
           '</claims></patent>')
    for c in independent_claims(xml):
        print("   [INDEPENDENT]", c)

    print("\n=== 3. resolver -> warhead SMARTS (all three formats) ===")
    for rep, kind in [("C=CC(=O)Nc1ccccc1", "smiles"),
                      ("InChI=1S/C9H9NO/c1-2-9(11)10-8-6-4-3-5-7-8/h2-7H,1H2,(H,10,11)", "inchi"),
                      ("N-(3-(5-chloro-1H-pyrrolo[2,3-b]pyridin-3-yl)phenyl)acrylamide", "iupac")]:
        s = resolve_to_smiles(rep, kind)
        print(f"   {kind:7} -> {str(s):46} warheads={has_warhead(s)}")
