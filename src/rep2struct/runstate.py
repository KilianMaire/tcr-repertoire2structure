from __future__ import annotations
import json
from pathlib import Path
from .schema import to_jsonable

class RunState:
    def __init__(self, run_dir):
        self.run_dir = Path(run_dir)
        self.run_dir.mkdir(parents=True, exist_ok=True)

    def path_for(self, name: str) -> Path:
        return self.run_dir / f"{name}.json"

    def stage_done(self, name: str) -> bool:
        return self.path_for(name).exists()

    def write_stage(self, name: str, obj) -> None:
        def enc(o):
            if isinstance(o, list):
                return [enc(x) for x in o]
            return to_jsonable(o)
        tmp = self.path_for(name).with_suffix(".json.tmp")
        tmp.write_text(json.dumps(enc(obj), indent=2))
        tmp.replace(self.path_for(name))

    def read_stage(self, name: str):
        return json.loads(self.path_for(name).read_text())
