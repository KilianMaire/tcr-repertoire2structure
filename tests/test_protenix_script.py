from rep2struct.tools import protenix_inputs, protenix_script

FASTA = (">A\nAAAA\n>B\nBBBB\n>C\nCCCC\n>D\nDDDD\n>E\nGILGFVFTL\n")


def _inputs():
    built = protenix_inputs.build(FASTA)
    return {f"c0_{k}": v for k, v in built.items()}


def test_script_is_bash_and_writes_inputs_and_folds():
    s = protenix_script.build(_inputs(), working_path="/scratch/run")
    assert s.startswith("#!/usr/bin/env bash")
    assert "set -euo pipefail" in s
    # both constructs are referenced
    assert "c0_cognate" in s and "c0_scramble" in s
    # the proven fold command and the local out/ layout QC reads
    assert "protenix pred -i inputs/" in s
    assert "-o out/" in s
    assert "--use_default_params true" in s
    # runs in the working path the user gave, no browser repatriation
    assert "/scratch/run" in s
    assert "files.download" not in s


def test_embeds_input_json_so_no_external_files_needed():
    s = protenix_script.build(_inputs())
    assert "GILGFVFTL" in s  # the peptide from the embedded record


def test_working_path_with_space_stays_intact_when_quoted():
    s = protenix_script.build(_inputs(), working_path="/scratch/my run")
    assert 'cd "/scratch/my run"' in s
