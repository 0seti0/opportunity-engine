"""OCSR — read DRAWN chemical structures off patent page images (the warheads that never appear in text;
the BI-HER2 case where Examples are 2D diagrams). Complements feed_build's text screen: when a patent reads
covalent (warhead mentions high) but no compound NAME resolves, its warheads are drawn -> come here.

Stack on this arm64 mac:
  - OCSR backend is pluggable (_ocsr_backend), default 'auto': prefers MolNexTR (Apache-2.0, ~93% recall on
    real USPTO patent images vs DECIMER ~41%) whenever it's importable AND a checkpoint is resolvable (the
    .venv-engine + .models/molnextr_best.pth setup), else DECIMER. OCSR_BACKEND=decimer is the kill switch.
    So .venv-engine runs use MolNexTR automatically; ephemeral `uv run` envs (no torch) use DECIMER. MolNexTR
    is the lead-resolution upgrade (the dossier's leads were 0/12 resolved on DECIMER). Graceful fallback.
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
_molnextr_model = None
_warned = False


def _model():
    global _predict
    if _predict is None:
        from DECIMER import predict_SMILES                    # heavy TF import — lazy, once
        _predict = predict_SMILES
    return _predict


def _molnextr_ckpt():
    """Resolve the MolNexTR checkpoint: OCSR_MOLNEXTR_CKPT if set+exists, else the .models/molnextr_best.pth
    symlink beside this module (the one-time-setup default). None if neither is present."""
    import os
    p = os.environ.get("OCSR_MOLNEXTR_CKPT")
    if p and os.path.exists(p):
        return p
    default = Path(__file__).parent / ".models" / "molnextr_best.pth"
    return str(default) if default.exists() else None


def _ocsr_backend():
    """Active OCSR backend. Default 'auto': prefer MolNexTR when it's importable AND a checkpoint resolves
    (_molnextr_ckpt), else DECIMER. OCSR_BACKEND=decimer forces DECIMER (kill switch). Best-available wins —
    .venv-engine picks MolNexTR; ephemeral runs (no torch) pick DECIMER. Pure logic, no model loaded."""
    import importlib.util
    import os
    if os.environ.get("OCSR_BACKEND", "auto").lower() == "decimer":
        return "decimer"
    if importlib.util.find_spec("MolNexTR") and _molnextr_ckpt():
        return "molnextr"
    return "decimer"


def _molnextr_smiles(png_path):
    """MolNexTR OCSR -> raw SMILES. Enable with:
        pip install git+https://github.com/CYF2000127/MolNexTR
        download molnextr_best.pth (HuggingFace CYF200127/MolNexTR, ~1.1GB), then
        export OCSR_MOLNEXTR_CKPT=/path/molnextr_best.pth OCSR_BACKEND=molnextr
    Checkpoint loaded once into the module-global model. API verified against the installed package: the class
    is MolNexTR.molnextr.molnextr(model_path, device); single-image call is .predict_final_results() returning
    a dict with 'predicted_smiles' (the README's top-level molnextr().prediction() is stale)."""
    global _molnextr_model
    if _molnextr_model is None:
        import torch
        from MolNexTR.molnextr import molnextr                # class lives in the same-named submodule
        _molnextr_model = molnextr(_molnextr_ckpt(), torch.device("cpu"))
    out = _molnextr_model.predict_final_results(str(png_path))
    return out.get("predicted_smiles") if isinstance(out, dict) else out


def _ocsr_raw(png_path):
    """Single-molecule PNG -> raw SMILES via the active backend, falling back to DECIMER on ANY MolNexTR
    error (missing install/checkpoint, load failure) so an OPS run never crashes on the optional dep."""
    global _warned
    if _ocsr_backend() == "molnextr":
        try:
            return _molnextr_smiles(png_path)
        except Exception as e:
            if not _warned:
                print(f"  [ocsr] MolNexTR unavailable ({type(e).__name__}: {e}); using DECIMER")
                _warned = True
    return _model()(str(png_path))


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
    raw = _ocsr_raw(png_path)
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
        raw = _ocsr_raw(d / "s.png")
        m = Chem.MolFromSmiles(raw) if raw else None
        smi = _clean(m) if m else None
        if smi:
            out.append(smi)
    return list(dict.fromkeys(out))


if __name__ == "__main__":
    import sys
    if "--selftest" in sys.argv:                          # backend selection only — no heavy model needed
        import importlib.util
        import os
        os.environ.pop("OCSR_BACKEND", None)
        avail = bool(importlib.util.find_spec("MolNexTR") and _molnextr_ckpt())
        assert _ocsr_backend() == ("molnextr" if avail else "decimer"), "auto must prefer MolNexTR iff available"
        os.environ["OCSR_BACKEND"] = "decimer"
        assert _ocsr_backend() == "decimer", "OCSR_BACKEND=decimer kill switch must force DECIMER"
        print(f"ocsr backend-selection self-check OK (auto -> {'molnextr' if avail else 'decimer'})")
        sys.exit()

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
