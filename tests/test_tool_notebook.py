import json
from rep2struct.tools.notebook import build_notebook


def test_notebook_is_valid_and_embeds_inputs():
    nb = build_notebook("af3", {"cognate": {"chains": {"E": "SII"}}})  # af3 is still unwired
    assert nb["nbformat"] == 4
    json.dumps(nb)  # serializable
    src = "".join(src_ for cell in nb["cells"] for src_ in cell["source"])
    assert "SII" in src                       # inputs embedded
    assert "TODO(live)" in src                 # live marker present


def test_notebook_scaffold_fails_loud_not_fake():
    nb = build_notebook("af3", {})             # af3 is not yet wired
    src = "".join(s for cell in nb["cells"] for s in cell["source"])
    assert "NotImplementedError" in src        # unwired notebook cannot fake a result


def test_tcrdock_notebook_is_wired_not_a_stub():
    inputs = {"c_cognate": {"row": {"organism": "human", "mhc_class": 1, "mhc": "A*02:01",
                                    "peptide": "GILGFVFTL", "va": "TRAV8-3*01", "ja": "TRAJ42*01",
                                    "cdr3a": "CAVGARGGSQGNLIF", "vb": "TRBV19*01",
                                    "jb": "TRBJ2-7*01", "cdr3b": "CASSTRAGVEQYF"}},
              "c_scramble": {"row": {"organism": "human", "mhc_class": 1, "mhc": "A*02:01",
                                     "peptide": "TFVFGLIGL", "va": "TRAV8-3*01", "ja": "TRAJ42*01",
                                     "cdr3a": "CAVGARGGSQGNLIF", "vb": "TRBV19*01",
                                     "jb": "TRBJ2-7*01", "cdr3b": "CASSTRAGVEQYF"}}}
    nb = build_notebook("tcrdock", inputs)
    json.dumps(nb)                                     # serializable
    src = "".join(s for cell in nb["cells"] for s in cell["source"])
    assert "live cell not yet validated" not in src    # not the fail-loud stub scaffold
    assert "TODO(live)" not in src
    assert "GILGFVFTL" in src and "TFVFGLIGL" in src   # both records embedded
    # the corrections that were required live to make TCRdock run:
    for marker in ("setup_for_alphafold", "run_prediction", "model_2_ptm_pae_1_2",
                   "v2.3.2", "jaxlib==0.3.25", "biopython==1.79", "cudatoolkit=11.1",
                   "LD_LIBRARY_PATH"):
        assert marker in src, f"missing recipe marker: {marker}"


def test_tcrdock_installs_jaxlib_before_af232_requirements():
    # Learned live 2026-07-09: jaxlib is NOT on PyPI (only Google's jax_releases index), and
    # chex inside the AF 2.3.2 requirements needs jaxlib. So jax+jaxlib (with the find-links
    # index) MUST be installed BEFORE `pip install -r af232.txt`, else the af232 resolve dies
    # with "No matching distribution found for jaxlib". Lock the order.
    inputs = {"k": {"row": {"organism": "human", "mhc_class": 1, "mhc": "A*02:01",
                            "peptide": "GILGFVFTL", "va": "TRAV8-3*01", "ja": "TRAJ42*01",
                            "cdr3a": "CAV", "vb": "TRBV19*01", "jb": "TRBJ2-7*01", "cdr3b": "CAS"}}}
    src = "".join(s for cell in build_notebook("tcrdock", inputs)["cells"] for s in cell["source"])
    assert src.index("jaxlib==0.3.25") < src.index("install -q -r /content/af232.txt")


def test_tcrdock_inverts_pae_for_verdict_binding():
    # interface PAE is LOWER = more recognized; verdict_binding treats HIGHER = more
    # recognized. The adapter must write score = -pae, else the verdict is inverted.
    inputs = {"k": {"row": {"organism": "human", "mhc_class": 1, "mhc": "A*02:01",
                            "peptide": "GILGFVFTL", "va": "TRAV8-3*01", "ja": "TRAJ42*01",
                            "cdr3a": "CAV", "vb": "TRBV19*01", "jb": "TRBJ2-7*01", "cdr3b": "CAS"}}}
    nb = build_notebook("tcrdock", inputs)
    src = "".join(s for cell in nb["cells"] for s in cell["source"])
    assert "score = -pae" in src


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


def test_affinetune_notebook_repatriates_the_scores():
    nb = build_notebook("affinetune", {"k": {"mhc": "M" * 200, "b2m": "B", "peptide": "GILGFVFTL"}})
    src = "".join(s for cell in nb["cells"] for s in cell["source"])
    assert "make_archive" in src and "affinetune_scores" in src
    assert "from google.colab import files" in src and "files.download" in src
    assert "DOWNLOAD_SKIPPED" in src


def test_tcrdock_notebook_repatriates_the_scores():
    inputs = {"k": {"row": {"organism": "human", "mhc_class": 1, "mhc": "A*02:01",
                            "peptide": "GILGFVFTL", "va": "TRAV8-3*01", "ja": "TRAJ42*01",
                            "cdr3a": "CAV", "vb": "TRBV19*01", "jb": "TRBJ2-7*01", "cdr3b": "CAS"}}}
    src = "".join(s for cell in build_notebook("tcrdock", inputs)["cells"] for s in cell["source"])
    assert "make_archive" in src and "tcrdock_scores" in src
    assert "files.download" in src and "DOWNLOAD_SKIPPED" in src


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


def test_protenix_notebook_is_wired_not_a_stub():
    inputs = {
        "c1_cognate": [{"name": "cognate", "sequences": [
            {"proteinChain": {"sequence": "GILGFVFTL", "count": 1, "id": ["E"]}}],
            "covalent_bonds": []}],
        "c1_scramble": [{"name": "scramble", "sequences": [
            {"proteinChain": {"sequence": "TFVFGLIGL", "count": 1, "id": ["E"]}}],
            "covalent_bonds": []}],
    }
    nb = build_notebook("protenix", inputs)
    json.dumps(nb)                                     # serializable
    src = "".join(s for cell in nb["cells"] for s in cell["source"])
    assert "live cell not yet validated" not in src    # not the fail-loud stub scaffold
    assert "TODO(live)" not in src
    assert "GILGFVFTL" in src and "TFVFGLIGL" in src   # both records embedded
    assert "c1_cognate" in src and "c1_scramble" in src
    for marker in ("pip install -q protenix", "protenix pred",
                   "protenix_base_default_v1.0.0", "--use_msa false"):
        assert marker in src, f"missing recipe marker: {marker}"


def test_protenix_notebook_repatriates_the_cifs():
    # The loop is only closed if the agent can pull the folded CIFs back. The notebook must
    # zip the out/ tree (preserving {cid}_cognate/{cid}_scramble in the paths) and hand it to
    # the browser via files.download so the Playwright executor captures it.
    nb = build_notebook("protenix", {"k": [{"name": "k", "sequences": [], "covalent_bonds": []}]})
    src = "".join(s for cell in nb["cells"] for s in cell["source"])
    assert "make_archive" in src and "protenix_folds" in src
    assert "from google.colab import files" in src and "files.download" in src
    # guarded so a non-Colab run does not die on the download
    assert "DOWNLOAD_SKIPPED" in src


def test_protenix_is_msa_free_matching_the_documented_reliable_run():
    # The run documented as reliable (docs/fold_qc_results.md) folded MSA-free after the
    # Protenix MSA server throttled; msa.py keeps the MSA out of the fold runtime. Lock it.
    nb = build_notebook("protenix", {"k": [{"name": "k", "sequences": [], "covalent_bonds": []}]})
    src = "".join(s for cell in nb["cells"] for s in cell["source"])
    assert "--use_msa false" in src and "--use_msa true" not in src


def test_inputs_cell_is_executable_python_with_bool_and_none():
    # JSON literals (false/null/true) are NOT valid Python; the INPUTS cell must
    # embed a Python literal so it executes in Jupyter/Colab without NameError.
    inputs = {"use_msa": False, "template": None, "pep": "SII"}
    nb = build_notebook("mhcfine", inputs)
    cell = next(c for c in nb["cells"] if "INPUTS =" in "".join(c["source"]))
    ns: dict = {}
    exec("".join(cell["source"]), ns)          # must not raise
    assert ns["INPUTS"] == inputs              # round-trips with False and None intact
