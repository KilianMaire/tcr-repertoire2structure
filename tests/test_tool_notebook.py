import json
from rep2struct.tools.notebook import build_notebook


def test_notebook_is_valid_and_embeds_inputs():
    nb = build_notebook("tcrdock", {"cognate": {"chains": {"E": "SII"}}})
    assert nb["nbformat"] == 4
    json.dumps(nb)  # serializable
    src = "".join(src_ for cell in nb["cells"] for src_ in cell["source"])
    assert "SII" in src                       # inputs embedded
    assert "TODO(live)" in src                 # live marker present


def test_notebook_scaffold_fails_loud_not_fake():
    nb = build_notebook("tcrdock", {})         # tcrdock is not yet wired
    src = "".join(s for cell in nb["cells"] for s in cell["source"])
    assert "NotImplementedError" in src        # unwired notebook cannot fake a result


def test_mhcfine_notebook_is_wired_not_a_stub():
    inputs = {"cognate": {"protein_sequence": "GSHSMRYFF", "peptide_sequence": "GILGFVFTL"},
              "scramble": {"protein_sequence": "GSHSMRYFF", "peptide_sequence": "TFVFGLIGL"}}
    nb = build_notebook("mhcfine", inputs)
    json.dumps(nb)                              # serializable
    src = "".join(s for cell in nb["cells"] for s in cell["source"])
    assert "NotImplementedError" not in src     # mhcfine is validated, not a stub
    assert "GILGFVFTL" in src and "TFVFGLIGL" in src   # both records embedded
    for marker in ("numpy<2", "kalign", "msa_run", "from src import", "inference"):
        assert marker in src, f"missing recipe marker: {marker}"


def test_mhcfine_deps_precede_import_so_no_restart():
    # The critical ordering learned live: numpy<2 + kalign must be installed BEFORE any
    # numpy/torch/src import, else np.string_ crashes and the kernel needs a restart.
    nb = build_notebook("mhcfine", {"cognate": {"protein_sequence": "M", "peptide_sequence": "SII"}})
    src = "".join(s for cell in nb["cells"] for s in cell["source"])
    assert src.index("numpy<2") < src.index("from src import")
    assert src.index("kalign") < src.index("import torch")


def test_inputs_cell_is_executable_python_with_bool_and_none():
    # JSON literals (false/null/true) are NOT valid Python; the INPUTS cell must
    # embed a Python literal so it executes in Jupyter/Colab without NameError.
    inputs = {"use_msa": False, "template": None, "pep": "SII"}
    nb = build_notebook("mhcfine", inputs)
    cell = next(c for c in nb["cells"] if "INPUTS =" in "".join(c["source"]))
    ns: dict = {}
    exec("".join(cell["source"]), ns)          # must not raise
    assert ns["INPUTS"] == inputs              # round-trips with False and None intact
