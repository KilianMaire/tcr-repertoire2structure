"""Select and render the Fig S5 gallery in one canonical orientation.

For each gallery epitope, pick the clonotype whose cognate is that epitope (a real
pairing, not a decoy fold), reconstructed (no poly-G stub), and most confident by
the TCR-to-peptide interface ipTM, fold it against its cognate, and render the whole
complex in the shared TCR-up frame (render_structure_pymol single_uniform). This
recovers and records the gallery provenance, which was not committed with the
original renders, and gives every panel the same orientation.

Runs in the PyMOL micromamba env:
  micromamba run -n pymol python scripts/render_gallery.py
"""
import glob
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from render_structure_pymol import render_single_uniform, render_peptide_inset

ROOT = Path(__file__).resolve().parent.parent
FIGS = ROOT / "paper/figures"

# panel letter, epitope, run, output stem (kept identical to the committed panels)
GALLERY = [
    ("a", "GILGFVFTL", "panel1", "_gallery1_GILGFVFTL.png"),
    ("b", "ELAGIGILTV", "panel1", "_gallery2_ELAGIGILTV.png"),
    ("c", "FLYALALLL", "panel1", "_gallery3_FLYALALLL.png"),
    ("d", "IVTDFSVIK", "hla_a1101", "_gallery4_IVTDFSVIK.png"),
]


def _is_stub(entry):
    return "GGGGGGGGGG" in entry.get("chain_b_seq", "")


def _iptm_tcrpep(summary):
    cpi = summary["chain_pair_iptm"]
    return max(cpi[0][4], cpi[1][4])


def _best_sample(pred_dir):
    """Return (cif_path, metric) for the sample with the best TCR-to-peptide ipTM."""
    best = None
    for js in sorted(glob.glob(str(pred_dir / "*summary_confidence_sample_*.json"))):
        metric = _iptm_tcrpep(json.load(open(js)))
        cif = js.replace("_summary_confidence_sample_", "_sample_").replace(".json", ".cif")
        if Path(cif).exists() and (best is None or metric > best[1]):
            best = (cif, metric)
    return best


def _select(epitope, run):
    """Best reconstructed clonotype whose cognate is `epitope`, by TCR-peptide ipTM."""
    manifest = json.load(open(ROOT.parent / f"runs/{run}/manifest.json"))
    ranked = []
    for cid, entry in manifest.items():
        if entry["cognate"] != epitope or _is_stub(entry):
            continue
        pred = ROOT.parent / f"runs/{run}/folds/{cid}__{epitope}/{cid}__{epitope}/seed_101/predictions"
        got = _best_sample(pred)
        if got:
            ranked.append((cid, got[0], got[1]))
    ranked.sort(key=lambda r: r[2], reverse=True)
    return ranked[0] if ranked else None


def main():
    print(f"{'panel':5} {'epitope':12} {'clonotype':14} {'iptm_TCRpep':>11}  file")
    for letter, epitope, run, out in GALLERY:
        sel = _select(epitope, run)
        if not sel:
            print(f"{letter:5} {epitope:12} {'NONE':14}")
            continue
        cid, cif, metric = sel
        render_single_uniform(cif, str(FIGS / out))
        render_peptide_inset(cif, str(FIGS / out.replace(".png", "_inset.png")))
        print(f"{letter:5} {epitope:12} {cid:14} {metric:11.3f}  {out}")


if __name__ == "__main__":
    main()
