"""OCSR — read DRAWN chemical structures off patent page images (the warheads that never appear in text;
the BI-HER2 case where Examples are 2D diagrams). Complements feed_build's text screen: when a patent reads
covalent (warhead mentions high) but no compound NAME resolves, its warheads are drawn -> come here.

Stack on this arm64 mac:
  - DECIMER (TensorFlow, models cached in ~/.data/DECIMER-V2) does SINGLE-MOLECULE recognition. Verified.
  - decimer-segmentation (the official page splitter) is INFEASIBLE here — it pins tensorflow<=2.15.1, which
    has no arm64 wheels. So we segment with a cv2 heuristic (dilate strokes into blobs, keep structure-sized,
    moderate-ink-density boxes) instead of the TF1 Mask-RCNN. Crude but dependency-free and arm64-clean.
  - macOS `sips` rasterizes the one-page OPS PDFs to PNG (no poppler/brew).
ponytail: cv2 heuristic segmentation, swap in a real segmenter (MolScribe / decimer-segmentation in a py3.10
or Docker env) if recall on dense scheme pages matters.
"""
import subprocess
import tempfile
from pathlib import Path

_predict = None


def _model():
    global _predict
    if _predict is None:
        from DECIMER import predict_SMILES                    # heavy TF import — lazy, once
        _predict = predict_SMILES
    return _predict


def pdf_to_png(pdf_bytes, dpi=300):
    """One-page PDF bytes -> PNG path (macOS sips). Lives in a temp dir the caller may leave to the OS."""
    d = Path(tempfile.mkdtemp())
    (d / "p.pdf").write_bytes(pdf_bytes)
    subprocess.run(["sips", "-s", "format", "png", str(d / "p.pdf"), "--out", str(d / "p.png")],
                   capture_output=True, check=True)
    return str(d / "p.png")


def image_to_smiles(png_path):
    """Single-molecule image -> canonical SMILES (or None)."""
    from rdkit import Chem
    raw = _model()(str(png_path))
    m = Chem.MolFromSmiles(raw) if raw else None
    return Chem.MolToSmiles(m) if m else None


def _crop_structures(png_path, max_crops=10):
    """Heuristic structure detection: dilate ink into blobs, keep structure-sized, line-drawing-density boxes.
    Thresholds tuned against real OPS patent pages — looser than first-guess (tight filters missed real
    structures: e.g. WO2026115265 p74 gave 0 tight vs 2 warheads loose)."""
    import cv2
    img = cv2.imread(png_path, cv2.IMREAD_GRAYSCALE)
    H, W = img.shape
    _, bw = cv2.threshold(img, 200, 255, cv2.THRESH_BINARY_INV)         # ink -> white
    merged = cv2.dilate(bw, cv2.getStructuringElement(cv2.MORPH_RECT, (20, 20)), iterations=1)
    boxes = []
    n, _, stats, _ = cv2.connectedComponentsWithStats(merged, 8)
    for i in range(1, n):
        x, y, w, h, _ = stats[i]
        if not (0.005 < (w * h) / (H * W) < 0.85 and 0.15 < w / h < 9):  # structure-sized, not a hairline/whole page
            continue
        if 0.015 < bw[y:y + h, x:x + w].mean() / 255 < 0.55:           # line drawing, not a solid text block
            boxes.append((w * h, img[max(0, y - 12):y + h + 12, max(0, x - 12):x + w + 12]))
    boxes.sort(key=lambda b: b[0], reverse=True)
    return [crop for _, crop in boxes[:max_crops]]


def _clean(mol):
    """Largest fragment as canonical SMILES, or None. Strips DECIMER noise fragments ('CC.CC...', atom salts)
    and real counter-ions alike; the warhead always sits on the main scaffold. Sub-drug-size => a noise crop."""
    from rdkit import Chem
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    big = max(frags, key=lambda f: f.GetNumHeavyAtoms()) if frags else mol
    return Chem.MolToSmiles(big) if big.GetNumHeavyAtoms() >= 8 else None


def page_structures(png_path):
    """Full patent page PNG -> list of canonical SMILES from the drawn structures on it. Crops the page (cv2),
    OCSRs each, keeps the largest fragment. Garbage survives twice-filtered: sub-drug-size crops are dropped
    here, and downstream the warhead SMARTS is itself strong (a random misread won't match an acrylamide)."""
    import cv2
    from rdkit import Chem
    out = []
    for crop in _crop_structures(png_path):
        d = Path(tempfile.mkdtemp())
        cv2.imwrite(str(d / "s.png"), crop)
        raw = _model()(str(d / "s.png"))
        m = Chem.MolFromSmiles(raw) if raw else None
        smi = _clean(m) if m else None
        if smi:
            out.append(smi)
    return list(dict.fromkeys(out))


if __name__ == "__main__":
    # round-trip self-check: render a known acrylamide, OCSR it back, confirm the warhead survives
    from rdkit import Chem
    from rdkit.Chem.Draw import MolToImage
    import patent_radar as pr
    d = Path(tempfile.mkdtemp())
    MolToImage(Chem.MolFromSmiles("C=CC(=O)Nc1ccccc1"), size=(400, 400)).save(str(d / "a.png"))
    smi = image_to_smiles(str(d / "a.png"))
    print("OCSR round-trip:", smi, "warheads:", pr.has_warhead(smi))
    assert smi and "acrylamide" in pr.has_warhead(smi), "DECIMER failed to recover the acrylamide warhead"
    print("OK")
