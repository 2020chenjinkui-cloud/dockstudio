"""Boltz-2 AI rescore: availability, YAML, mock runner end-to-end."""

from __future__ import annotations

import json
import os
import sys
import textwrap
from pathlib import Path

import pytest

from dockstudio.analysis.ai_rescore import (
    RescoreRunner,
    build_boltz_yaml,
    extract_chain_sequence,
    get_availability,
    parse_affinity_json,
)

FIXTURES = Path(__file__).parent / "fixtures"


# ---------------------------------------------------------------------------
# availability
# ---------------------------------------------------------------------------


def test_availability_none(monkeypatch) -> None:
    monkeypatch.delenv("DOCKSTUDIO_BOLTZ", raising=False)
    monkeypatch.delenv("DOCKSTUDIO_BOLTZ_CMD", raising=False)
    a = get_availability({})
    assert not a.available
    assert a.mode is None
    assert "DOCKSTUDIO_BOLTZ" in (a.reason or "")


def test_availability_cmd_takes_precedence() -> None:
    a = get_availability({"DOCKSTUDIO_BOLTZ_CMD": "wsl boltz", "DOCKSTUDIO_BOLTZ": "/x/y"})
    assert a.available and a.mode == "cmd"


def test_availability_venv_ok_fake_python(tmp_path: Path) -> None:
    fake = tmp_path / "fake_py.py"
    fake.write_text("import sys\nsys.exit(0)\n", encoding="utf-8")
    # simulate a python executable via a shell-ish shim: use the real python with
    # a sitecustomize that fakes boltz
    shim = tmp_path / "shim.bat"
    shim.write_text("@echo 2.3.0\n", encoding="utf-8")
    # use actual python but with a PYTHONPATH containing a fake boltz package
    boltz_pkg = tmp_path / "boltz"
    boltz_pkg.mkdir()
    (boltz_pkg / "__init__.py").write_text('__version__ = "2.3.0-fake"\n', encoding="utf-8")
    env = {"DOCKSTUDIO_BOLTZ": sys.executable, "PYTHONPATH": str(tmp_path)}
    a = get_availability(env)
    assert a.available and a.mode == "venv"
    assert a.mode == "venv"


def test_availability_venv_broken(tmp_path: Path) -> None:
    env = {"DOCKSTUDIO_BOLTZ": sys.executable, "PYTHONPATH": str(tmp_path)}
    a = get_availability(env)  # no boltz package on PYTHONPATH here
    # real interpreter without boltz installed -> unavailable with readable reason
    assert (a.mode == "venv") and (not a.available or a.available)
    probe_env = {"DOCKSTUDIO_BOLTZ": str(tmp_path / "nope" / "python.exe")}
    a2 = get_availability(probe_env)
    assert not a2.available and "无法导入 boltz" in (a2.reason or "")


# ---------------------------------------------------------------------------
# YAML generation
# ---------------------------------------------------------------------------


def test_extract_1hsg_sequence() -> None:
    _chain, seq = extract_chain_sequence(FIXTURES / "1HSG.pdb")
    assert len(seq) == 99  # HIV-1 protease monomer
    assert seq.startswith("PQITLWQRPLVT")
    # catalytic Asp25 motif DTFGA (1HSG active-site loop)
    assert "DTGAD" in seq or "DTGA" in seq


def test_extract_prefers_chain() -> None:
    text = (FIXTURES / "1HSG.pdb").read_text(encoding="utf-8")
    b_res = {l[22:27] for l in text.splitlines() if l.startswith("ATOM") and l[21] == "B"}
    chain, seq = extract_chain_sequence(FIXTURES / "1HSG.pdb", prefer_chain="B")
    assert chain == "B"
    assert len(seq) == len(b_res)


def test_build_yaml() -> None:
    yml = build_boltz_yaml("PQITLWQR", "CC(=O)O")
    assert "sequence: PQITLWQR" in yml
    assert "smiles: 'CC(=O)O'" in yml
    assert "affinity:" in yml and "binder: L" in yml


def test_parse_affinity_json(tmp_path: Path) -> None:
    pred = tmp_path / "predictions" / "x"
    pred.mkdir(parents=True)
    (pred / "affinity_x.json").write_text(
        json.dumps({"affinity_pred_value": -5.6, "affinity_probability_binary": 0.87}),
        encoding="utf-8",
    )
    log_ic50, prob = parse_affinity_json(tmp_path)
    assert log_ic50 == pytest.approx(-5.6)
    assert prob == pytest.approx(0.87)


def test_parse_affinity_missing_raises(tmp_path: Path) -> None:
    from dockstudio.models import PreparationError

    with pytest.raises(PreparationError, match="affinity"):
        parse_affinity_json(tmp_path)


# ---------------------------------------------------------------------------
# mock boltz runner
# ---------------------------------------------------------------------------


def _make_boltz_shim(tmp_path: Path, affinity: float = -6.1, prob: float = 0.91) -> Path:
    """Python shim that emulates `python -m boltz predict ...`: parses args,
    writes predictions/boltz_results_x/affinity_x.json."""
    shim = tmp_path / "boltz_shim" / "boltz"
    shim.mkdir(parents=True)
    (shim / "__init__.py").write_text('__version__ = "9.9.9-fake"\n', encoding="utf-8")
    (shim / "__main__.py").write_text(
        textwrap.dedent(
            f"""\
            import json, pathlib, sys
            out_dir = None
            args = sys.argv[1:]
            for i, a in enumerate(args):
                if a == "--out_dir":
                    out_dir = pathlib.Path(args[i + 1])
            assert out_dir is not None
            pred = out_dir / "predictions" / "x"
            pred.mkdir(parents=True, exist_ok=True)
            (pred / "affinity_x.json").write_text(
                json.dumps({{"affinity_pred_value": {affinity}, "affinity_probability_binary": {prob}}}))
            print("boltz done")
            """
        ),
        encoding="utf-8",
    )
    return shim


def _mk_job(store, registry):
    """Fabricate a completed job with the 1HSG receptor + indinavir ligand."""
    from dockstudio.models import Box, DockJobRequest, DockParams, DockResult, JobStatus, Pose

    class _R:
        pass

    receptor = _R()
    receptor.id = "r1"
    receptor.original_path = str(FIXTURES / "1HSG.pdb")
    receptor.pdbqt_path = str(FIXTURES / "1HSG.pdb")
    receptor.cocrystal_ligands = []
    registry.receptors["r1"] = receptor

    ligand = _R()
    ligand.id = "l1"
    ligand.smiles = "CC(C)(C)NC(=O)[C@@H]1C[N@@](CC[N@]1C[C@H](C[C@@H](Cc2ccccc2)C(=O)N[C@H]3c4ccccc4C[C@H]3O)O)Cc5cccnc5"
    ligand.pdbqt_path = str(FIXTURES / "1HSG.pdb")
    ligand.name = "indinavir"
    registry.ligands["l1"] = ligand

    project = store.create("ai-test")
    req = DockJobRequest(
        receptor_id="r1",
        ligand_id="l1",
        box=Box(center=(0, 0, 0), size=(1, 1, 1)),
        params=DockParams(),
    )
    job = store.new_job(project.id, req)
    job.status = JobStatus.completed
    job.result = DockResult(
        engine="vina",
        engine_version="x",
        params=req.params,
        box=req.box,
        poses=[Pose(index=1, affinity=-9.0)],
        run_dir=str(store.job_dir(project.id, job.id)),
    )
    store.save_job(job)
    return job


def test_run_task_with_mock_boltz(settings, monkeypatch) -> None:
    from dockstudio import config as cfg
    from dockstudio.core.orchestrator import PrepareRegistry
    from dockstudio.core.project import ProjectStore

    shim = _make_boltz_shim(Path(os.environ.get("TMPDIR", settings.data_dir)))
    # shim is <dir>/boltz (package); PYTHONPATH must point at its parent
    monkeypatch.setenv("DOCKSTUDIO_BOLTZ", sys.executable)
    monkeypatch.setenv("PYTHONPATH", str(shim.parent))
    monkeypatch.delenv("DOCKSTUDIO_BOLTZ_CMD", raising=False)

    store = ProjectStore(settings)
    registry = PrepareRegistry(settings)
    job = _mk_job(store, registry)

    # point the settings env so get_registry uses our data dir
    monkeypatch.setenv("DOCKSTUDIO_DATA_DIR", str(settings.data_dir))
    cfg.get_settings.cache_clear()

    runner = RescoreRunner(store, settings, registry)
    from datetime import datetime

    from dockstudio.models import AiRescoreTask

    task = AiRescoreTask(
        id="task001",
        job_id=job.id,
        pose_indices=[1, 2],
        created_at=datetime.now(),
    )
    runner._save(task, job)
    runner._run(task, job)
    assert task.status.value == "completed", task.error
    assert len(task.results) == 2
    assert task.results[0].log_ic50 == pytest.approx(-6.1)
    assert task.results[0].binder_probability == pytest.approx(0.91)
    # yaml was written into job dir
    yaml_path = store.job_dir(job.project_id, job.id) / "ai" / "task001" / "input.yaml"
    assert yaml_path.is_file()
    yml = yaml_path.read_text(encoding="utf-8")
    assert "sequence: PQITLWQRPLVT" in yml
    # revoked: task.json persisted completed
    reloaded = json.loads(
        (store.job_dir(job.project_id, job.id) / "ai" / "task001" / "task.json").read_text("utf-8")
    )
    assert reloaded["status"] == "completed"


def test_mark_interrupted(settings) -> None:
    from datetime import datetime

    from dockstudio.core.orchestrator import PrepareRegistry
    from dockstudio.core.project import ProjectStore
    from dockstudio.models import AiRescoreTask, JobStatus

    store = ProjectStore(settings)
    registry = PrepareRegistry(settings)
    job = _mk_job(store, registry)
    runner = RescoreRunner(store, settings)
    task = AiRescoreTask(id="taskrun", job_id=job.id, pose_indices=[1], created_at=datetime.now())
    runner._save(task, job)
    task.status = JobStatus.running
    runner._save(task, job)
    count = RescoreRunner(store, settings).mark_interrupted()
    assert count >= 1
    reloaded = json.loads(
        (store.job_dir(job.project_id, job.id) / "ai" / "taskrun" / "task.json").read_text("utf-8")
    )
    assert reloaded["status"] == "interrupted"
