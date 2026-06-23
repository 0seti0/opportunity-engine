"""Step 0+1 of the covalent early-warning radar — the Day-1 patent screen (no SureChEMBL wait).

  independent_claims(xml) -> the independent claims (where the CORE invention lives; we search only here)
  warhead_in_claims(text) -> syntax-bound warhead regex + exclusion list  (catches NAMED warheads)
  resolve_to_smiles(rep)  -> a compound (SMILES/InChI | IUPAC name | image) -> SMILES, preferring the most
                             reliable representation (RDKit | OPSIN | DECIMER)
  has_warhead(smiles)     -> covalent-warhead SMARTS                       (catches DRAWN warheads)

Runs on raw patent XML at publication. SureChEMBL becomes a later cross-check, not the source.
"""
import re
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

from rdkit import Chem

# --- covalent-warhead SMARTS (acyclic constraints exclude steroid-enone / fumarate false positives) -------
WARHEAD_SMARTS = {
    "acrylamide":       "[CX3;!R]=[CX3;!R]C(=O)[NX3]",
    "chloroacetamide":  "[Cl]C[CX3](=O)[NX3]",
    "propiolamide":     "[CX2]#[CX2]C(=O)[NX3]",
    "vinylsulfonamide": "[CX3;!R]=[CX3;!R]S(=O)(=O)[NX3]",
}
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
    """ST.36-style patent XML: a <claim> is DEPENDENT if it references another (<claim-ref> or 'claim N')."""
    root = ET.fromstring(xml_text)
    out = []
    for cl in root.iter():
        if cl.tag.split('}')[-1] == "claim":
            txt = "".join(cl.itertext())
            dependent = any(c.tag.split('}')[-1] == "claim-ref" for c in cl.iter()) \
                or re.search(r"\bclaim\s+\d+\b", txt, re.I)
            if not dependent:
                out.append(re.sub(r"\s+", " ", txt).strip())
    return out


def _opsin(name):
    """IUPAC name -> SMILES. Production: bundled OPSIN jar (needs a JRE). Here: hosted fallback so it runs."""
    try:
        from py2opsin import py2opsin                          # local OPSIN (no network) — preferred in prod
        s = py2opsin(name)
        if s:
            return s
    except Exception:
        pass
    for url in (f"https://opsin.ch/opsin/{urllib.parse.quote(name)}.smi",
                f"https://cactus.nci.nih.gov/chemical/structure/{urllib.parse.quote(name)}/smiles"):
        try:
            s = urllib.request.urlopen(url, timeout=20).read().decode().strip()
            if s and Chem.MolFromSmiles(s):
                return s
        except Exception:
            continue
    return None


def resolve_to_smiles(rep, kind):
    """rep -> canonical SMILES. kind: smiles | inchi | iupac | image. Prefer the most reliable available."""
    if kind == "smiles":
        m = Chem.MolFromSmiles(rep)
    elif kind == "inchi":
        m = Chem.MolFromInchi(rep)
    elif kind == "iupac":
        s = _opsin(rep)
        m = Chem.MolFromSmiles(s) if s else None
    else:                                                       # image -> DECIMER OCSR (ML fallback) — prod stub
        return None
    return Chem.MolToSmiles(m) if m else None


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
