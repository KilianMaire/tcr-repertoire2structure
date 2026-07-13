from __future__ import annotations

import json

TOOL = "protenix"

# The proven Protenix invocation, identical to tools/notebook._protenix_notebook cell 3.
_FOLD_CMD = ("protenix pred -i inputs/{key}.json -o out/{key} -s 101 "
             "-n protenix_base_default_v1.0.0 --use_default_params true")


def build(inputs: dict, working_path: str = ".") -> str:
    """A self contained bash script that folds each embedded Protenix record on a machine
    the user has a shell on. INPUTS is {key: <protenix prediction JSON>} (cognate + scramble,
    keys prefixed by clonotype id). CIFs land under out/{key} in working_path, the same
    {cid}_cognate / {cid}_scramble layout the beta V-domain to peptide QC calibrates on. MSA free by
    design (mirrors the documented reliable Protenix path); no browser repatriation because
    the outputs are already local."""
    embedded = json.dumps(inputs)
    lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        f'cd "{working_path}"',
        "pip install -q protenix",
        "mkdir -p inputs out",
        "python - <<'PY'",
        "import json, os",
        f"INPUTS = json.loads({embedded!r})",
        "os.makedirs('inputs', exist_ok=True)",
        "for key, obj in INPUTS.items():",
        "    json.dump(obj, open(f'inputs/{key}.json', 'w'))",
        "print('wrote', len(INPUTS), 'inputs:', sorted(INPUTS))",
        "PY",
        "for f in inputs/*.json; do",
        '  key=$(basename "$f" .json)',
        f"  {_FOLD_CMD.format(key='$key')}",
        "done",
        "echo DONE",
    ]
    return "\n".join(lines) + "\n"
