"""Ray-traced cartoon (ribbon) renders of a predicted TCR-pMHC complex.

Runs inside the PyMOL micromamba env, not the project .venv:
  micromamba run -n pymol python scripts/render_structure_pymol.py <cif> <out_dir>

Produces two ray-traced PNGs (whole complex, and looking down the groove) with the
five chains coloured. The compositor scripts/plot_structure_3d.py assembles them
into the annotated figure. Colours match that figure's legend.
"""
import sys
from pymol import cmd

CIF = sys.argv[1] if len(sys.argv) > 1 else "runs/panel1/folds/b97fd808da3f__GILGFVFTL/b97fd808da3f__GILGFVFTL/seed_101/predictions/b97fd808da3f__GILGFVFTL_sample_1.cif"
OUT = sys.argv[2] if len(sys.argv) > 2 else "paper/figures"

COLORS = {"A": "0x0072B2", "B": "0x009E73", "C": "0xE69F00", "D": "0xCC79A7", "E": "0xD55E00"}

cmd.reinitialize()
cmd.load(CIF, "m")
cmd.hide("everything")
cmd.dss("m")                          # assign secondary structure for cartoon
cmd.show("cartoon", "m")
cmd.set("cartoon_transparency", 0.0)
cmd.set("cartoon_fancy_helices", 1)
for ch, col in COLORS.items():
    cmd.color(col, f"m and chain {ch}")
# peptide antigen: emphasise as sticks over its cartoon
cmd.show("sticks", "m and chain E")
cmd.set("stick_radius", 0.35, "m and chain E")
cmd.set("cartoon_transparency", 0.0, "m and chain E")

cmd.bg_color("white")
cmd.set("ray_opaque_background", 1)
cmd.set("antialias", 2)
cmd.set("ambient", 0.35)
cmd.set("specular", 0.2)
cmd.set("ray_shadows", 0)
cmd.orient("m")

cmd.ray(1500, 1500)
cmd.png(f"{OUT}/_struct_view1.png", dpi=200)

cmd.turn("x", 90)                     # look down onto the peptide groove
cmd.ray(1500, 1500)
cmd.png(f"{OUT}/_struct_view2.png", dpi=200)
print("rendered", f"{OUT}/_struct_view1.png", f"{OUT}/_struct_view2.png")
