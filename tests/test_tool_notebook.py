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
    for marker in ("np.string_", "np.sum", "kalign", "msa_run", "from src import", "inference"):
        assert marker in src, f"missing recipe marker: {marker}"


def test_mhcfine_keeps_stock_numpy2_no_downgrade():
    # Learned live 2026-07-09: DOWNGRADING numpy (numpy<2) poisons numpy's own compiled
    # mtrand.so on today's Colab image (dtype ABI wall). The robust recipe KEEPS stock
    # numpy 2 and shims the removed AF2-era aliases instead. Lock the lesson: no downgrade.
    nb = build_notebook("mhcfine", {"cognate": {"protein_sequence": "M", "peptide_sequence": "SII"}})
    src = "".join(s for cell in nb["cells"] for s in cell["source"])
    assert "numpy<2" not in src and "numpy==1" not in src


def test_mhcfine_shim_precedes_import_so_no_restart():
    # The critical ordering: the numpy-2 compat shim (np.string_, np.sum-on-generator, ...)
    # and kalign must be in place BEFORE `from src import`, else the AF2-derived code hits
    # the removed np.string_ / np.sum(generator) and the run dies mid-fold.
    nb = build_notebook("mhcfine", {"cognate": {"protein_sequence": "M", "peptide_sequence": "SII"}})
    src = "".join(s for cell in nb["cells"] for s in cell["source"])
    assert src.index("np.string_") < src.index("from src import")
    assert src.index("kalign") < src.index("from src import")


def test_affinetune_notebook_is_wired_not_a_stub():
    inputs = {"c_cognate": {"mhc": "GSHSMRYFF" * 25, "b2m": "IQRTPKIQV", "peptide": "GILGFVFTL"},
              "c_scramble": {"mhc": "GSHSMRYFF" * 25, "b2m": "IQRTPKIQV", "peptide": "TFVFGLIGL"}}
    nb = build_notebook("affinetune", inputs)
    json.dumps(nb)                                     # serializable
    src = "".join(s for cell in nb["cells"] for s in cell["source"])
    assert "live cell not yet validated" not in src    # not the fail-loud stub scaffold
    assert "TODO(live)" not in src
    assert "GILGFVFTL" in src and "TFVFGLIGL" in src   # both records embedded
    for marker in ("run_prediction", "model_2_ptm_ft_pae", "--ignore_identities",
                   "cudatoolkit=11.1", "LD_LIBRARY_PATH", "[:175]"):
        assert marker in src, f"missing recipe marker: {marker}"


def test_affinetune_inverts_pae_for_verdict_binding():
    # run_prediction pae is LOWER = presented; verdict_binding treats HIGHER = presented.
    # The adapter must write score = -pae, else the verdict is inverted. Lock the direction.
    nb = build_notebook("affinetune", {"k": {"mhc": "M" * 200, "b2m": "B", "peptide": "GILGFVFTL"}})
    src = "".join(s for cell in nb["cells"] for s in cell["source"])
    assert "score = -pae" in src


def test_affinetune_ninemer_only_is_fail_loud():
    # Only 9-mers are calibrated; a non-9-mer must raise, not silently fold against a
    # length-mismatched template. The guard lives in the embedded loop.
    nb = build_notebook("affinetune", {"k": {"mhc": "M" * 200, "b2m": "B", "peptide": "SIINFEKL"}})
    src = "".join(s for cell in nb["cells"] for s in cell["source"])
    assert "NotImplementedError" in src and "only 9-mers wired" in src


def test_inputs_cell_is_executable_python_with_bool_and_none():
    # JSON literals (false/null/true) are NOT valid Python; the INPUTS cell must
    # embed a Python literal so it executes in Jupyter/Colab without NameError.
    inputs = {"use_msa": False, "template": None, "pep": "SII"}
    nb = build_notebook("mhcfine", inputs)
    cell = next(c for c in nb["cells"] if "INPUTS =" in "".join(c["source"]))
    ns: dict = {}
    exec("".join(cell["source"]), ns)          # must not raise
    assert ns["INPUTS"] == inputs              # round-trips with False and None intact
