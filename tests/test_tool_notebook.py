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
    nb = build_notebook("mhcfine", {})
    src = "".join(s for cell in nb["cells"] for s in cell["source"])
    assert "NotImplementedError" in src        # unwired notebook cannot fake a result


def test_inputs_cell_is_executable_python_with_bool_and_none():
    # JSON literals (false/null/true) are NOT valid Python; the INPUTS cell must
    # embed a Python literal so it executes in Jupyter/Colab without NameError.
    inputs = {"use_msa": False, "template": None, "pep": "SII"}
    nb = build_notebook("mhcfine", inputs)
    cell = next(c for c in nb["cells"] if "INPUTS =" in "".join(c["source"]))
    ns: dict = {}
    exec("".join(cell["source"]), ns)          # must not raise
    assert ns["INPUTS"] == inputs              # round-trips with False and None intact
