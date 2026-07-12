"""Ray-traced cartoon (ribbon) renders of a predicted TCR-pMHC complex.

Runs inside the PyMOL micromamba env, not the project .venv:
  micromamba run -n pymol python scripts/render_structure_pymol.py complex     <cif> <out_dir>
  micromamba run -n pymol python scripts/render_structure_pymol.py groove      <cif> <out.png>
  micromamba run -n pymol python scripts/render_structure_pymol.py groove_conf <cif> <out.png>
  micromamba run -n pymol python scripts/render_structure_pymol.py interface   <cif> <out.png>

Chain order from the pipeline: A=TCR alpha, B=TCR beta, C=MHC heavy, D=beta-2
microglobulin, E=peptide. Colours match figstyle.CHAIN and the Fig 1 legend.

Renders are deliberately clean (no baked-in text). Text annotation is added by the
matplotlib compositors, so labels can be placed and kept uniform across figures.

  complex     -> _struct_view1.png (TCR-up side view) + _struct_view2.png (groove top-down)
  groove      -> one top-down groove crop (MHC helices amber + peptide sticks)
  groove_conf -> groove crop, peptide coloured by per-residue pLDDT, residues labelled
  interface   -> the recognition interface: TCR loops that contact the peptide, as sticks
"""
import math
import sys
from pymol import cmd

COLORS = {"A": "0x0072B2", "B": "0x009E73", "C": "0xE69F00", "D": "0xCC79A7", "E": "0xD55E00"}
ONE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
       "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
       "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
       "TYR": "Y", "VAL": "V"}


def _base_scene(cif):
    cmd.reinitialize()
    cmd.load(cif, "m")
    cmd.hide("everything")
    cmd.dss("m")
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)
    cmd.set("antialias", 2)
    cmd.set("ambient", 0.35)
    cmd.set("specular", 0.2)
    cmd.set("ray_shadows", 0)
    cmd.set("cartoon_fancy_helices", 1)


def _color_all():
    for ch, col in COLORS.items():
        cmd.color(col, f"m and chain {ch}")


def _sub(a, b):
    return [a[i] - b[i] for i in range(3)]


def _norm(v):
    n = math.sqrt(sum(c * c for c in v)) or 1.0
    return [c / n for c in v]


def _cross(a, b):
    return [a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]]


def _dot(a, b):
    return sum(a[i] * b[i] for i in range(3))


def _orient_tcr_up():
    """Canonical side view built from geometry, not from get_view. Up = the docking
    axis (peptide centroid to TCR centroid); right = the peptide N-to-C axis made
    orthogonal to up; the camera looks along their cross product. Sets the rotation
    directly so the TCR variable domains sit above the peptide-MHC platform."""
    cmd.orient("m")
    tcr = cmd.centerofmass("m and (chain A or chain B)")
    pep = cmd.centerofmass("m and chain E")
    cas = cmd.get_model("m and chain E and name CA").atom
    n_term, c_term = cas[0].coord, cas[-1].coord

    up = _norm(_sub(tcr, pep))                        # TCR is above the peptide
    right0 = _norm(_sub(c_term, n_term))              # peptide runs left to right
    right = _norm(_sub(right0, [up[i] * _dot(right0, up) for i in range(3)]))
    view_dir = _norm(_cross(right, up))               # camera z, toward viewer

    v = cmd.get_view()
    new = list(right) + list(up) + list(view_dir) + list(v[9:])
    cmd.set_view(new)


def _orient_groove_topdown():
    """Deterministic top-down groove view, built from geometry so every structure
    lands in the same reading frame (unlike cmd.orient, which picks a per-structure
    camera). Image-right = the peptide long axis (first principal component of the CA
    trace), signed N-to-C so the N terminus is always left and C right; toward-viewer
    = the peptide-plane normal (smallest principal component), signed to the TCR side
    so the bulge points up out of the page. The peptide therefore lies flat and reads
    left-to-right the same way in every panel."""
    # PyMOL picks the natural, well-spread top-down of the peptide (long axis in
    # plane, thinnest dimension toward the viewer); we only fix its 4-fold sign
    # ambiguity so the frame is identical across structures.
    cmd.orient("m and chain E")
    tcr = cmd.centerofmass("m and (chain A or chain B)")
    cas = cmd.get_model("m and chain E and name CA").atom
    n_term, c_term = cas[0].coord, cas[-1].coord

    def _screen_x(p):                                 # image-x of a model point (column 0)
        v = cmd.get_view()
        return v[0] * p[0] + v[3] * p[1] + v[6] * p[2]

    def _screen_z(p):                                 # toward-viewer of a model point (column 2)
        v = cmd.get_view()
        return v[2] * p[0] + v[5] * p[1] + v[8] * p[2]

    if _screen_x(c_term) < _screen_x(n_term):         # want C to the right, N to the left
        cmd.turn("y", 180)
    if _screen_z(tcr) < _screen_z(cmd.centerofmass("m and chain E")):  # view from the TCR side
        cmd.turn("x", 180)


def _label_termini():
    """Prominent N and C tags just beyond the peptide termini so the reading
    direction is explicit, not left to the (occludable) residue numbers. The tags
    are pushed outward along the N-to-C axis so they clear the sticks and the
    per-residue labels."""
    cas = cmd.get_model("m and chain E and name CA").atom
    n_term, c_term = cas[0].coord, cas[-1].coord
    axis = _norm(_sub(c_term, n_term))
    off = [4.0 * a for a in axis]
    ends = (("N", n_term, [-off[i] for i in range(3)]),
            ("C", c_term, off))
    for tag, base, delta in ends:
        name = f"term_{tag}"
        cmd.pseudoatom(name, pos=[base[i] + delta[i] for i in range(3)], label=tag)
        cmd.hide("nonbonded", name)
        cmd.set("label_size", 36, name)
        cmd.set("label_color", "0x111111", name)


def _set_camera(right, up, view_dir):
    """Set the camera rotation from an explicit orthonormal basis. PyMOL stores the
    first nine view elements column-major, so the basis vectors go in as columns."""
    v = cmd.get_view()
    m = [right[0], up[0], view_dir[0],
         right[1], up[1], view_dir[1],
         right[2], up[2], view_dir[2]]
    cmd.set_view(list(m) + list(v[9:]))


def _orient_complex_uniform():
    """Canonical whole-complex view shared by every gallery panel: the TCR sits above
    the peptide-MHC platform (image-up = the docking axis), the peptide runs left to
    right (image-right = the peptide N-to-C axis), and the camera looks from the side.
    Built with the column-major set_view convention so the basis is applied, not its
    transpose, which is what let earlier renders tilt differently per structure."""
    cmd.orient("m")
    tcr = cmd.centerofmass("m and (chain A or chain B)")
    pep = cmd.centerofmass("m and chain E")
    cas = cmd.get_model("m and chain E and name CA").atom
    n_term, c_term = cas[0].coord, cas[-1].coord

    up = _norm(_sub(tcr, pep))                        # TCR above the platform
    right0 = _norm(_sub(c_term, n_term))              # peptide N to C, made horizontal
    right = _norm(_sub(right0, [up[i] * _dot(right0, up) for i in range(3)]))
    view_dir = _norm(_cross(right, up))
    _set_camera(right, up, view_dir)


def render_single_uniform(cif, out_png):
    """Whole-complex gallery render in the canonical TCR-up frame."""
    _base_scene(cif)
    cmd.show("cartoon", "m")
    _color_all()
    cmd.show("sticks", "m and chain E")
    cmd.set("stick_radius", 0.35, "m and chain E")
    _orient_complex_uniform()
    cmd.zoom("m", 3)
    cmd.ray(1450, 1500)
    cmd.png(out_png, dpi=200)
    print("rendered single_uniform", out_png)


def render_complex(cif, out_dir):
    _base_scene(cif)
    cmd.show("cartoon", "m")
    _color_all()
    cmd.show("sticks", "m and chain E")
    cmd.set("stick_radius", 0.35, "m and chain E")
    _orient_tcr_up()
    cmd.zoom("m", 3)
    cmd.ray(1450, 1600)
    cmd.png(f"{out_dir}/_struct_view1.png", dpi=200)

    # view 2: look straight down onto the peptide groove, TCR and b2m removed.
    cmd.hide("everything", "m and (chain A or chain B or chain D)")
    cmd.orient("(m and chain C within 12 of (m and chain E)) or (m and chain E)")
    cmd.zoom("(m and chain C within 12 of (m and chain E)) or (m and chain E)", 4)
    cmd.ray(1500, 1150)
    cmd.png(f"{out_dir}/_struct_view2.png", dpi=200)
    print("rendered complex", out_dir)


def render_groove(cif, out_png):
    """Top-down groove: MHC heavy chain cartoon (amber) + peptide sticks (orange)."""
    _base_scene(cif)
    cmd.show("cartoon", "m and chain C")
    cmd.color(COLORS["C"], "m and chain C")
    cmd.set("cartoon_transparency", 0.35, "m and chain C")
    cmd.show("sticks", "m and chain E")
    cmd.color(COLORS["E"], "m and chain E")
    cmd.set("stick_radius", 0.32, "m and chain E")
    cmd.orient("(m and chain C within 12 of (m and chain E)) or (m and chain E)")
    cmd.zoom("(m and chain C within 10 of (m and chain E)) or (m and chain E)", 5)
    cmd.ray(1200, 950)
    cmd.png(out_png, dpi=200)
    print("rendered groove", out_png)


def render_peptide_inset(cif, out_png):
    """Tight top-down zoom on the peptide for a gallery inset: the MHC groove is a
    recessive translucent grey so the orange peptide sticks are the clear focus."""
    _base_scene(cif)
    cmd.show("cartoon", "m and chain C")
    cmd.color("grey80", "m and chain C")
    cmd.set("cartoon_transparency", 0.65, "m and chain C")
    cmd.show("sticks", "m and chain E")
    cmd.color(COLORS["E"], "m and chain E")
    cmd.set("stick_radius", 0.4, "m and chain E")
    cmd.orient("(m and chain C within 10 of (m and chain E)) or (m and chain E)")
    cmd.zoom("m and chain E", 3.5)
    cmd.ray(1150, 900)
    cmd.png(out_png, dpi=200)
    print("rendered peptide_inset", out_png)


def render_groove_conf(cif, out_png, pmin=50, pmax=95):
    """Top-down groove, peptide sticks coloured by per-residue pLDDT (B-factor),
    each residue labelled one-letter + position. High pLDDT -> blue, low -> red."""
    _base_scene(cif)
    cmd.show("cartoon", "m and chain C")
    cmd.color("grey80", "m and chain C")
    cmd.set("cartoon_transparency", 0.6, "m and chain C")
    cmd.show("sticks", "m and chain E")
    cmd.set("stick_radius", 0.35, "m and chain E")
    cmd.spectrum("b", "red_white_blue", "m and chain E", minimum=pmin, maximum=pmax)
    cmd.set("label_size", 19)
    cmd.set("label_color", "black")
    cmd.set("float_labels", 1)
    for at in cmd.get_model("m and chain E and name CA").atom:
        aa = ONE.get(at.resn, "X")
        cmd.label(f"m and chain E and name CA and resi {at.resi}", f'"{aa}{at.resi}"')
    _orient_groove_topdown()
    _label_termini()
    cmd.zoom("m and chain E", 4.5)
    cmd.ray(1300, 950)
    cmd.png(out_png, dpi=200)
    print("rendered groove_conf", out_png)


def render_interface(cif, out_png):
    """The recognition interface: TCR loops within 5A of the peptide shown as sticks
    over a translucent TCR cartoon, with the peptide sticks and a ghost MHC groove."""
    _base_scene(cif)
    cmd.show("cartoon", "m and (chain A or chain B)")
    cmd.color(COLORS["A"], "m and chain A")
    cmd.color(COLORS["B"], "m and chain B")
    cmd.set("cartoon_transparency", 0.55, "m and (chain A or chain B)")
    cmd.show("cartoon", "m and chain C")
    cmd.color("grey80", "m and chain C")
    cmd.set("cartoon_transparency", 0.75, "m and chain C")
    cmd.show("sticks", "m and chain E")
    cmd.color(COLORS["E"], "m and chain E")
    cmd.set("stick_radius", 0.32, "m and chain E")
    cmd.select("contacts", "m and (chain A or chain B) within 5 of (m and chain E)")
    cmd.show("sticks", "contacts")
    cmd.set("stick_radius", 0.26, "contacts")
    cmd.orient("(m and chain E) or contacts")
    cmd.zoom("(m and chain E) or contacts", 5)
    cmd.ray(1300, 1050)
    cmd.png(out_png, dpi=200)
    print("rendered interface", out_png)


def render_single(cif, out_png):
    """Just the TCR-up whole-complex view (canonical orientation) to a named PNG,
    for the supplementary gallery of multiple complexes."""
    _base_scene(cif)
    cmd.show("cartoon", "m")
    _color_all()
    cmd.show("sticks", "m and chain E")
    cmd.set("stick_radius", 0.35, "m and chain E")
    _orient_tcr_up()
    cmd.zoom("m", 3)
    cmd.ray(1450, 1500)
    cmd.png(out_png, dpi=200)
    print("rendered single", out_png)


DISPATCH = {"complex": lambda a: render_complex(a[0], a[1]),
            "groove": lambda a: render_groove(a[0], a[1]),
            "groove_conf": lambda a: render_groove_conf(a[0], a[1]),
            "interface": lambda a: render_interface(a[0], a[1]),
            "single": lambda a: render_single(a[0], a[1]),
            "single_uniform": lambda a: render_single_uniform(a[0], a[1]),
            "peptide_inset": lambda a: render_peptide_inset(a[0], a[1])}

if __name__ == "__main__":
    mode = sys.argv[1]
    if mode not in DISPATCH:
        raise SystemExit(f"unknown mode {mode!r}; use one of {sorted(DISPATCH)}")
    DISPATCH[mode](sys.argv[2:])
