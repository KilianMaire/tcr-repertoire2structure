"""Assemble the manuscript into a single Word document with every figure embedded.

Combines paper/main.md (main text plus the six main figures), the graphical
abstract, paper/methods.md, and paper/supplementary/supplementary.md (with the five
supplementary figures inlined), then calls pandoc to write paper/Repertoire2Structure.docx.
Images are the committed PNGs, so no GPU or PyMOL is needed.

Usage: python scripts/build_docx.py   (needs pandoc on PATH)
"""
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PAPER = ROOT / "paper"
FIGDIR = "paper/figures"          # path as pandoc will resolve it (run from ROOT)
OUT = PAPER / "Repertoire2Structure.docx"
WIDTH = "{width=6.3in}"

FIG_FILE = {
    "1": "fig1_structure.png", "2": "fig2_validation.png", "3": "fig3_retrieval.png",
    "4": "fig4_confidence_variance.png", "5": "fig5_mhc_presentation.png",
    "6": "fig6_two_axis_map.png",
}
SUPP_FILE = {
    "S1": "figS1_retrieval_strata.png", "S2": "figS2_groove_confidence.png",
    "S3": "figS3_chain_pair_iptm.png", "S4": "figS4_reproducibility.png",
    "S5": "figS5_complex_gallery.png",
}


def _img(fname):
    return f"![]({FIGDIR}/{fname}){WIDTH}"


def build_main():
    text = (PAPER / "main.md").read_text()
    body, figures = text.split("\n## Figures", 1)

    # graphical abstract just before the Introduction
    ga = (f"## Graphical abstract\n\n{_img('graphical_abstract.png')}\n\n"
          "**Graphical abstract.** Structural confidence reads MHC-peptide "
          "presentation (cognate versus scramble groove, AUROC up to 0.99) but not "
          "TCR-peptide recognition (held-out retrieval 0.61, 11 of 18, p 0.24).\n\n")
    body = body.replace("## Introduction", ga + "## Introduction", 1)

    # rebuild the Figures section: image above each caption, as a normal paragraph
    out = ["## Figures\n"]
    for num, cap in re.findall(r"- \*\*Figure (\d)\.\*\* (.*?)(?=\n- \*\*Figure|\nSee |\Z)",
                               figures, flags=re.S):
        cap = " ".join(cap.split())
        out.append(f"\n{_img(FIG_FILE[num])}\n\n**Figure {num}.** {cap}\n")
    return body + "\n".join(out)


def build_supp():
    text = (PAPER / "supplementary/supplementary.md").read_text()
    # inline each supplementary figure right after its header
    def repl(m):
        tag = m.group(1)
        return f"{m.group(0)}\n\n{_img(SUPP_FILE[tag])}\n"
    return re.sub(r"### Figure (S\d)\..*", repl, text)


def main():
    methods = (PAPER / "methods.md").read_text()
    combined = "\n\n".join([
        build_main(),
        "\\newpage\n\n" + methods,
        "\\newpage\n\n" + build_supp(),
    ])
    md_path = PAPER / "_manuscript_combined.md"
    md_path.write_text(combined)

    cmd = ["pandoc", str(md_path.relative_to(ROOT)), "-o", str(OUT.relative_to(ROOT)),
           "--resource-path", str(ROOT), "--from", "gfm+tex_math_dollars",
           "--metadata", "title=Structural confidence in TCR-pMHC prediction reads "
           "presentation, not recognition"]
    subprocess.run(cmd, cwd=ROOT, check=True)
    md_path.unlink()
    print("wrote", OUT.relative_to(ROOT))


if __name__ == "__main__":
    main()
