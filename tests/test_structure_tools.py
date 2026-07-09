from rep2struct import structure_tools as st


def test_default_is_protenix():
    d = st.get_default()
    assert d.name == "protenix" and d.is_default


def test_registry_has_the_five_v1_tools():
    names = {t.name for t in st.REGISTRY}
    assert names == {"protenix", "af3", "mhcfine", "tcrdock", "affinetune"}


def test_affinetune_is_binding_score_the_rest_structure():
    by = {t.name: t for t in st.REGISTRY}
    assert by["affinetune"].output_type == "binding_score"
    assert by["protenix"].output_type == "structure"


def test_tcrdock_output_type_stays_structure_despite_binding_qc():
    # tcrdock's qc_metric became binding_score (single-chain output, scored by interface PAE),
    # but its output_type MUST stay "structure": it is the TCR:pMHC docking tool and has to
    # remain selectable for structure jobs. Flipping output_type would drop it from
    # tools_for(output_needed="structure") and wrongly add it to binding_score coverage.
    by = {t.name: t for t in st.REGISTRY}
    assert by["tcrdock"].output_type == "structure"
    assert "tcrdock" in {t.name for t in st.tools_for(1, True, "human", "structure")}
    assert "tcrdock" not in {t.name for t in st.tools_for(1, True, "human", "binding_score")}


def test_tools_for_class_ii_structure_returns_protenix_not_mhcfine():
    got = {t.name for t in st.tools_for(2, has_tcr=True, species="human", output_needed="structure")}
    assert "protenix" in got and "mhcfine" not in got  # mhcfine is class I only


def test_tools_for_binding_score_class_ii_returns_affinetune():
    got = {t.name for t in st.tools_for(2, has_tcr=False, species="mouse", output_needed="binding_score")}
    assert got == {"affinetune"}


def test_is_covered_true_for_structure_false_when_no_match():
    assert st.is_covered(1, True, "human", "structure") is True
    assert st.is_covered(1, True, "human", "binding_score") is False  # has TCR + wants presentation score -> affinetune (needs_tcr=False) does not match


def test_as_dicts_is_json_safe():
    import json
    json.dumps(st.as_dicts())  # must not raise


def test_qc_metric_per_tool():
    by = {t.name: t.qc_metric for t in st.REGISTRY}
    # tcrdock folds a structure but its single-chain output is scored by the peptide<->TCR
    # interface PAE via verdict_binding, so its qc_metric is binding_score (its output_type
    # stays "structure" so it remains selectable for structure jobs; see the registry note).
    assert by == {"protenix": "cdr3_peptide", "tcrdock": "binding_score",
                  "mhcfine": "peptide_groove", "affinetune": "binding_score",
                  "af3": "cdr3_peptide"}


def test_qc_metric_for_defaults():
    assert st.qc_metric_for("mhcfine") == "peptide_groove"
    assert st.qc_metric_for("unknown") == "cdr3_peptide"


def test_as_dicts_exposes_qc_metric():
    assert all("qc_metric" in d for d in st.as_dicts())
