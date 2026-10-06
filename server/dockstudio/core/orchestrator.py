"""Run orchestration + prepared-entity registry for the M0 synchronous flow."""

from __future__ import annotations

import json
import shutil
import uuid
from datetime import datetime
from pathlib import Path

from ..config import Settings, get_settings
from ..logging import get_logger
from ..models import DockResult, DockStudioError, PreparedLigand, PreparedReceptor

log = get_logger("core.orchestrator")


class PrepareRegistry:
    """JSON registry mapping prepared ids to their model payloads."""

    def __init__(self, settings: Settings):
        self.settings = settings
        self.path = settings.prepared_dir / "prepare_registry.json"
        self.receptors: dict[str, PreparedReceptor] = {}
        self.ligands: dict[str, PreparedLigand] = {}
        self.load()

    @property
    def receptors_dir(self) -> Path:
        d = self.settings.prepared_dir / "receptors"
        d.mkdir(parents=True, exist_ok=True)
        return d

    @property
    def ligands_dir(self) -> Path:
        d = self.settings.prepared_dir / "ligands"
        d.mkdir(parents=True, exist_ok=True)
        return d

    def load(self) -> None:
        if not self.path.is_file():
            return
        try:
            data = json.loads(self.path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            log.warning("could not load prepare registry %s: %s", self.path, exc)
            return
        for payload in data.get("receptors", []):
            receptor = PreparedReceptor(**payload)
            self.receptors[receptor.id] = receptor
        for payload in data.get("ligands", []):
            ligand = PreparedLigand(**payload)
            self.ligands[ligand.id] = ligand

    def save(self) -> None:
        payload = {
            "receptors": [r.model_dump() for r in self.receptors.values()],
            "ligands": [lig.model_dump() for lig in self.ligands.values()],
        }
        tmp = self.path.with_suffix(".tmp")
        tmp.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
        tmp.replace(self.path)

    def add_receptor(self, receptor: PreparedReceptor) -> None:
        self.receptors[receptor.id] = receptor
        self.save()

    def add_ligand(self, ligand: PreparedLigand) -> None:
        self.ligands[ligand.id] = ligand
        self.save()


def new_run_dir(settings: Settings) -> Path:
    """Create runs/<timestamp>_<8-hex>/ for one docking job."""
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")  # noqa: DTZ005  # local time for humans
    for _ in range(5):
        run_dir = settings.runs_dir / f"{stamp}_{uuid.uuid4().hex[:8]}"
        try:
            run_dir.mkdir(parents=True)
            return run_dir
        except FileExistsError:
            continue
    raise DockStudioError("could not allocate a run directory")


def stage_inputs(result: DockResult, receptor_pdbqt: Path, ligand_pdbqt: Path) -> DockResult:
    """Copy prepared inputs into the run dir so runs are self-contained."""
    run_dir = Path(result.run_dir)
    shutil.copy2(receptor_pdbqt, run_dir / "receptor.pdbqt")
    shutil.copy2(ligand_pdbqt, run_dir / "ligand.pdbqt")
    return result


def get_registry(settings: Settings | None = None) -> PrepareRegistry:
    return PrepareRegistry(settings or get_settings())
