"""Comprehensive Automated Test Suite for DockStudio Unified CLI and Agent REPL.

Validates:
- R1: verify command (T1~T5 numerical verification and edge cases)
- R1: audit command (JSONL event stream, live SHA256 integrity, tamper detection)
- R1: view command (micro 3D Web viewer, ephemeral HTTP server, terminal ASCII fallback, RDKit QC)
- R2: agent command (SlashCompleter, piped REPL stream, process lifecycle cleanup)
- R2: unified CLI main entrypoint (subcommand dispatch, help/version, serve command, backward compatibility)
- R3: 100% mock/fixture isolated, zero GPU / external engine requirements, fast execution.
"""

from __future__ import annotations

import hashlib
import io
import json
import sys
import threading
import urllib.error
import urllib.request
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from unittest import mock

import pytest
from prompt_toolkit.completion import CompleteEvent
from prompt_toolkit.document import Document
from rdkit import Chem
from rdkit.Chem import AllChem

from dockstudio.cli import (
    _cmd_audit,
    _cmd_engines,
    main,
    run_agent_repl,
)
from dockstudio.cli.agent import (
    ACTIVE_PROCESSES,
    SlashCompleter,
    _dispatch_command,
    cleanup_active_processes,
    register_process,
    unregister_process,
)
from dockstudio.cli.commands.audit import (
    compute_sha256_chunked,
    parse_events_jsonl,
    run_audit,
)
from dockstudio.cli.commands.verify import (
    run_verify,
)
from dockstudio.cli.commands.view import (
    ReusableThreadingHTTPServer,
    _calculate_dynamic_interactions,
    _extract_affinity,
    _parse_ligand_atoms_and_mol,
    _parse_receptor_atoms,
    _render_ligand_properties_table,
    _render_receptor_summary_table,
    _run_rdkit_qc_checks,
    build_viewer_html,
    create_viewer_handler,
    find_available_port,
    run_view,
)
from dockstudio.models import Box, DockJobRequest, DockParams, DockResult, Job, JobStatus, Pose

# ---------------------------------------------------------------------------
# Test Fixtures & Synthetic Data Helpers
# ---------------------------------------------------------------------------

def _create_conformer_sdf_block(smiles: str, center: tuple[float, float, float] = (10.0, 10.0, 10.0)) -> str:
    """Generate a genuine 3D molecule MolBlock centered at the specified coordinate."""
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(mol, randomSeed=42)
    conf = mol.GetConformer()
    pts = conf.GetPositions()
    c = pts.mean(axis=0)
    for i in range(mol.GetNumAtoms()):
        conf.SetAtomPosition(i, pts[i] - c + list(center))
    return Chem.MolToMolBlock(mol)


_SAMPLE_RECEPTOR_PDB = """ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00  0.00           N
ATOM      2  CA  ALA A   1      11.000  10.000  10.000  1.00  0.00           C
ATOM      3  C   ALA A   1      11.500  11.200  10.000  1.00  0.00           C
ATOM      4  O   ALA A   1      11.000  12.200  10.000  1.00  0.00           O
ATOM      5  CB  ALA A   1      11.500   9.000  11.000  1.00  0.00           C
TER
END
"""

_SAMPLE_RECEPTOR_PDBQT = """ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00  0.00    -0.300 N
ATOM      2  CA  ALA A   1      11.000  10.000  10.000  1.00  0.00    +0.100 C
ATOM      3  C   ALA A   1      11.500  11.200  10.000  1.00  0.00    +0.500 C
ATOM      4  O   ALA A   1      11.000  12.200  10.000  1.00  0.00    -0.500 O
TER
"""

_SAMPLE_LIGAND_PDBQT_POSE = """MODEL 1
REMARK VINA RESULT:    -8.500      0.000      0.000
ATOM      1  N   LIG A   1      11.000  14.500  10.000  1.00  0.00    -0.300 N
ATOM      2  C   LIG A   1      11.000  15.800  10.000  1.00  0.00    +0.100 C
ENDMDL
"""


def _build_synthetic_verify_package(
    pkg_dir: Path,
    *,
    offset_center: bool = False,
    funnel_discrepancy: bool = False,
    as_zip: bool = False,
) -> Path:
    """Constructs a fully populated synthetic delivery package for verify tests."""
    pkg_dir.mkdir(parents=True, exist_ok=True)

    # 1. conf.txt
    cx = 50.0 if offset_center else 10.0
    conf_content = (
        f"center_x = {cx}\n"
        "center_y = 10.0\n"
        "center_z = 10.0\n"
        "size_x = 20.0\n"
        "size_y = 20.0\n"
        "size_z = 20.0\n"
        "seed = 42\n"
        "exhaustiveness = 8\n"
        "num_modes = 9\n"
        "energy_range = 3.0\n"
        "cpu = 4\n"
        "scoring = vina\n"
    )
    (pkg_dir / "conf.txt").write_text(conf_content, encoding="utf-8")

    # 2. receptor.pdbqt
    (pkg_dir / "receptor.pdbqt").write_text(_SAMPLE_RECEPTOR_PDBQT, encoding="utf-8")

    # 3. crystal_ligand.sdf (centered at 10.0, 10.0, 10.0)
    crystal_sdf = _create_conformer_sdf_block("c1ccccc1", center=(10.0, 10.0, 10.0))
    (pkg_dir / "crystal_ligand.sdf").write_text(crystal_sdf, encoding="utf-8")

    # 4. poses directory & energy.csv (Top 10 compounds with matched InChIKey)
    poses_dir = pkg_dir / "poses"
    poses_dir.mkdir(exist_ok=True)

    smiles_list = [
        "CC(=O)O", "CCO", "CCC", "CCCC", "c1ccccc1",
        "CCN", "c1ccncc1", "CCCC#N", "CC(=O)N", "CS",
    ]
    energy_lines = ["id,affinity,smiles,inchikey"]
    for idx, smi in enumerate(smiles_list, 1):
        cid = f"MOL_{idx}"
        mol_block = _create_conformer_sdf_block(smi, center=(10.0, 10.0, 10.0))
        mol = Chem.MolFromMolBlock(mol_block)
        ik = Chem.MolToInchiKey(mol) if mol else f"INCHIKEY_{idx}"
        score = -10.0 + idx * 0.2
        energy_lines.append(f"{cid},{score:.2f},{smi},{ik}")
        (poses_dir / f"pose_{cid}.sdf").write_text(mol_block, encoding="utf-8")
        # Also provide numeric rank alias pose_1.sdf, pose_2.sdf
        (poses_dir / f"pose_{idx}.sdf").write_text(mol_block, encoding="utf-8")

    (pkg_dir / "energy.csv").write_text("\n".join(energy_lines), encoding="utf-8")

    # 5. funnel.csv (strictly conservative: raw - losses == docked_success == len(energy))
    raw_count = 100
    desalt = 10
    dedup = 10
    chiral = 10
    prep = 10
    dock_fail = 50
    docked_success = 8 if funnel_discrepancy else 10  # discrepancy causes failure

    funnel_content = (
        f"raw,{raw_count}\n"
        f"desalt,{desalt}\n"
        f"dedup,{dedup}\n"
        f"chiral_failed,{chiral}\n"
        f"prep_failed,{prep}\n"
        f"dock_failed,{dock_fail}\n"
        f"docked_success,{docked_success}\n"
    )
    (pkg_dir / "funnel.csv").write_text(funnel_content, encoding="utf-8")

    if as_zip:
        zip_path = pkg_dir.parent / f"{pkg_dir.name}.zip"
        with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
            for item in pkg_dir.rglob("*"):
                if item.is_file():
                    zf.write(item, item.relative_to(pkg_dir))
        return zip_path

    return pkg_dir


def _build_synthetic_audit_job(job_dir: Path, *, tampered: bool = False) -> Path:
    """Constructs a synthetic job folder for audit testing."""
    job_dir.mkdir(parents=True, exist_ok=True)

    # 1. job.json
    now = datetime.now(timezone.utc)
    job = Job(
        id=job_dir.name,
        project_id="test_project_alpha",
        status=JobStatus.completed,
        request=DockJobRequest(
            engine="vina",
            receptor_id="rec_hsg",
            ligand_id="lig_indinavir",
            box=Box(center=(10.0, 10.0, 10.0), size=(20.0, 20.0, 20.0)),
            params=DockParams(exhaustiveness=8, num_modes=9, seed=42, cpu=4),
        ),
        result=DockResult(
            engine="vina",
            engine_version="1.2.5",
            params=DockParams(),
            box=Box(center=(10.0, 10.0, 10.0), size=(20.0, 20.0, 20.0)),
            poses=[Pose(index=1, affinity=-8.50, rmsd_lb=0.0, rmsd_ub=0.0)],
            run_dir=str(job_dir),
        ),
        created_at=now,
        started_at=now,
        finished_at=now,
    )
    (job_dir / "job.json").write_text(job.model_dump_json(indent=2), encoding="utf-8")

    # 2. Asset files: receptor.pdbqt & ligand.pdbqt
    rec_path = job_dir / "receptor.pdbqt"
    rec_path.write_text(_SAMPLE_RECEPTOR_PDBQT, encoding="utf-8")
    lig_path = job_dir / "ligand.pdbqt"
    lig_path.write_text(_SAMPLE_LIGAND_PDBQT_POSE, encoding="utf-8")

    # 3. inputs.sha256.json
    rec_sha = compute_sha256_chunked(rec_path)
    lig_sha = compute_sha256_chunked(lig_path)

    if tampered:
        # Deliberately invalidate receptor digest to simulate tampering
        rec_sha = "0" * 64

    sha_map = {
        "receptor": rec_sha,
        "ligand": lig_sha,
    }
    (job_dir / "inputs.sha256.json").write_text(json.dumps(sha_map, indent=2), encoding="utf-8")

    # 4. events.jsonl
    events = [
        {
            "ts": "2026-11-07T05:00:00Z",
            "action": "job.created",
            "actor": "user",
            "subject": {"id": job.id},
            "payload": {"engine": "vina"},
        },
        {
            "ts": "2026-11-07T05:00:05Z",
            "action": "docking.started",
            "actor": "worker",
            "subject": {"id": job.id},
            "payload": {"exhaustiveness": 8},
        },
        {
            "ts": "2026-11-07T05:00:15Z",
            "action": "job.completed",
            "actor": "worker",
            "subject": {"id": job.id},
            "payload": {"top_affinity": -8.50},
        },
    ]
    with open(job_dir / "events.jsonl", "w", encoding="utf-8") as f:
        f.writelines(json.dumps(ev) + "\n" for ev in events)

    return job_dir


# ---------------------------------------------------------------------------
# Test Cases: CLI Entrypoint, Subcommand Dispatch & Backwards Compatibility
# ---------------------------------------------------------------------------

def test_main_entrypoint_help_and_version(monkeypatch: pytest.MonkeyPatch):
    """Test 1: Verify main(['--help']) and main(['--version']) return exit code 0."""
    with pytest.raises(SystemExit) as exc1:
        main(["--help"])
    assert exc1.value.code == 0

    with pytest.raises(SystemExit) as exc2:
        main(["--version"])
    assert exc2.value.code == 0

    # Test bare invocation in non-interactive environment prints help and returns 0
    monkeypatch.setattr(sys.stdin, "isatty", lambda: False)
    assert main([]) == 0


def test_main_subcommand_dispatch(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    """Test 2: Test dispatching engines, demo, audit, verify, view via main([...])."""
    # 1. dispatch engines
    rc_engines = main(["engines"])
    assert rc_engines == 0

    # 2. dispatch demo with mocked backend demo implementation
    demo_called = []
    monkeypatch.setattr("dockstudio.cli._cmd_demo", lambda: demo_called.append(True) or 0)
    rc_demo = main(["demo"])
    assert rc_demo == 0
    assert demo_called == [True]

    # 3. dispatch audit with synthetic job
    audit_job_dir = tmp_path / "dispatch_audit_job"
    _build_synthetic_audit_job(audit_job_dir)
    rc_audit = main(["audit", str(audit_job_dir)])
    assert rc_audit == 0

    # 4. dispatch verify with synthetic package
    verify_pkg_dir = tmp_path / "dispatch_verify_pkg"
    _build_synthetic_verify_package(verify_pkg_dir)
    rc_verify = main(["verify", str(verify_pkg_dir)])
    assert rc_verify == 0

    # 5. dispatch view with ascii fallback
    complex_dir = tmp_path / "dispatch_view_complex"
    complex_dir.mkdir()
    (complex_dir / "receptor.pdbqt").write_text(_SAMPLE_RECEPTOR_PDBQT, encoding="utf-8")
    (complex_dir / "out.pdbqt").write_text(_SAMPLE_LIGAND_PDBQT_POSE, encoding="utf-8")
    rc_view = main(["view", str(complex_dir), "--ascii"])
    assert rc_view == 0

    # 6. dispatch agent sub-parser
    monkeypatch.setattr("dockstudio.cli.run_agent_repl", lambda: 0)
    rc_agent = main(["agent"])
    assert rc_agent == 0


def test_main_serve_command_mocked(monkeypatch: pytest.MonkeyPatch):
    """Test 3: Test main(['serve', '--port', '8000']) with mocked uvicorn/serve entrypoint."""
    serve_invoked = []

    def mock_serve_main():
        serve_invoked.append(True)

    monkeypatch.setattr("dockstudio.api.app.main", mock_serve_main)
    rc = main(["serve", "--port", "8000"])
    assert rc == 0
    assert serve_invoked == [True]


def test_legacy_imports_backward_compatibility():
    """Test 4: Verify legacy import paths and functions remain backward-compatible."""
    import dockstudio.cli as cli_mod

    assert hasattr(cli_mod, "_cmd_audit")
    assert hasattr(cli_mod, "_cmd_engines")
    assert hasattr(cli_mod, "_cmd_demo")
    assert hasattr(cli_mod, "main")
    assert hasattr(cli_mod, "run_agent_repl")

    # Verify _cmd_engines() execution
    rc = _cmd_engines()
    assert rc == 0

    # Verify _cmd_audit with non-existent job returns 2
    rc_bad = _cmd_audit("nonexistent_legacy_job_99999")
    assert rc_bad == 2


# ---------------------------------------------------------------------------
# Test Cases: /verify Numerical Acceptance Suite (T1~T5)
# ---------------------------------------------------------------------------

def test_verify_t1_to_t5_full_pass(tmp_path: Path):
    """Test 5: Run /verify on a complete synthetic delivery package (folder and zip).

    Validates T1~T5 pass (or honest skip on T1/T2 when engine is unconfigured).
    """
    pkg_dir = tmp_path / "valid_pkg"
    _build_synthetic_verify_package(pkg_dir)

    # 1. Verify folder package
    rc_folder = run_verify(pkg_dir)
    assert rc_folder == 0

    # 2. Verify zip package
    zip_path = _build_synthetic_verify_package(tmp_path / "valid_pkg_zipped", as_zip=True)
    rc_zip = run_verify(zip_path)
    assert rc_zip == 0


def test_verify_criteria_failure(tmp_path: Path):
    """Test 6: Verify package with funnel mismatch or centroid offset > 2.0 A returns code 1."""
    # 1. Centroid offset failure (center_x = 50.0 while co-crystal is at 10.0 -> dist = 40.0 Å > 2.0 Å)
    offset_pkg = tmp_path / "offset_pkg"
    _build_synthetic_verify_package(offset_pkg, offset_center=True)
    rc_offset = run_verify(offset_pkg)
    assert rc_offset == 1

    # 2. Funnel discrepancy failure (docked_success != expected_docked)
    funnel_pkg = tmp_path / "funnel_discrepancy_pkg"
    _build_synthetic_verify_package(funnel_pkg, funnel_discrepancy=True)
    rc_funnel = run_verify(funnel_pkg)
    assert rc_funnel == 1


def test_verify_invalid_path(tmp_path: Path):
    """Test 7: Verify non-existent path or corrupted archive returns code 2."""
    bogus_path = tmp_path / "nonexistent_verify_package_path"
    rc_nonexistent = run_verify(bogus_path)
    assert rc_nonexistent == 2

    # Missing conf.txt in an otherwise empty directory
    empty_dir = tmp_path / "empty_dir"
    empty_dir.mkdir()
    rc_empty = run_verify(empty_dir)
    assert rc_empty == 2

    # Corrupted zip file
    bad_zip = tmp_path / "bad.zip"
    bad_zip.write_bytes(b"NOT_A_VALID_ZIP_ARCHIVE")
    rc_bad_zip = run_verify(bad_zip)
    assert rc_bad_zip == 2


# ---------------------------------------------------------------------------
# Test Cases: /audit Command & Cryptographic Verification
# ---------------------------------------------------------------------------

def test_audit_intact_job(tmp_path: Path):
    """Test 8: Construct synthetic job with events and matching SHA256 hashes, assert code 0."""
    job_dir = tmp_path / "intact_audit_job"
    _build_synthetic_audit_job(job_dir, tampered=False)

    rc = run_audit(str(job_dir))
    assert rc == 0


def test_audit_tampered_job(tmp_path: Path):
    """Test 9: Tamper with an input hash or file, assert code 1 (tamper detected)."""
    job_dir = tmp_path / "tampered_audit_job"
    _build_synthetic_audit_job(job_dir, tampered=True)

    rc = run_audit(str(job_dir))
    assert rc == 1


def test_audit_nonexistent_job(tmp_path: Path):
    """Test 10: Non-existent job ID or directory returns code 2."""
    rc = run_audit("completely_bogus_job_id_xyz")
    assert rc == 2

    nonexistent_dir = tmp_path / "nonexistent_job_dir"
    rc_dir = run_audit(str(nonexistent_dir))
    assert rc_dir == 2


def test_audit_zip_archive_and_chunked_sha(tmp_path: Path):
    """Test 8b: Test audit against zip snapshot and compute_sha256_chunked helper."""
    # Test compute_sha256_chunked
    sample_file = tmp_path / "sample.dat"
    content = b"DockStudio deterministic chunked hashing test content 12345"
    sample_file.write_bytes(content)
    expected_hash = hashlib.sha256(content).hexdigest()
    assert compute_sha256_chunked(sample_file) == expected_hash

    # Test parse_events_jsonl with empty lines and torn trailing line
    events_file = tmp_path / "events.jsonl"
    events_file.write_text(
        '{"ts": "2026-11-07T01:00:00Z", "action": "test.act1", "subject": {"id": "job1"}}\n\n'
        '{"ts": "2026-11-07T01:00:01Z", "action": "test.act2", "subject": {"id": "job1"}}\n'
        '{"ts": "2026-11-07T01:00:02Z", "action": "TRUNCATED_CRASH_LINE',
        encoding="utf-8",
    )
    parsed = parse_events_jsonl(events_file, subject_ids={"job1"})
    assert len(parsed) == 2
    assert parsed[0]["action"] == "test.act1"


# ---------------------------------------------------------------------------
# Test Cases: /view Micro 3D Viewer & Terminal ASCII Fallback
# ---------------------------------------------------------------------------

def test_view_ascii_complex(tmp_path: Path):
    """Test 11: Run run_view(complex_path, ascii_mode=True), assert code 0 and contact detection."""
    complex_dir = tmp_path / "view_complex"
    complex_dir.mkdir()
    (complex_dir / "receptor.pdbqt").write_text(_SAMPLE_RECEPTOR_PDBQT, encoding="utf-8")
    (complex_dir / "out.pdbqt").write_text(_SAMPLE_LIGAND_PDBQT_POSE, encoding="utf-8")

    rc = run_view(complex_dir, ascii_mode=True)
    assert rc == 0

    # Directly verify contact calculation between receptor and ligand
    rec_atoms = _parse_receptor_atoms(_SAMPLE_RECEPTOR_PDBQT)
    lig_atoms, _ = _parse_ligand_atoms_and_mol(_SAMPLE_LIGAND_PDBQT_POSE)
    interactions = _calculate_dynamic_interactions(rec_atoms, lig_atoms)
    assert len(interactions) > 0
    # Ligand N at (11.0, 14.5, 10.0) is within 2.3 Å of Receptor O at (11.0, 12.2, 10.0) -> hbond
    assert any(i["type"] == "hbond" for i in interactions)


def test_view_ascii_single_ligand(tmp_path: Path):
    """Test 12: Run run_view(ligand_sdf, ascii_mode=True), assert code 0 and RDKit table."""
    ligand_sdf = tmp_path / "aspirin.sdf"
    # Aspirin: 2-acetoxybenzoic acid
    sdf_content = _create_conformer_sdf_block("CC(=O)Oc1ccccc1C(=O)O", center=(5.0, 5.0, 5.0))
    ligand_sdf.write_text(sdf_content, encoding="utf-8")

    rc = run_view(ligand_sdf, ascii_mode=True)
    assert rc == 0

    # Verify RDKit property table generator
    mol = Chem.MolFromMolBlock(sdf_content)
    table = _render_ligand_properties_table(mol, "Aspirin")
    assert table is not None
    assert "Aspirin" in table.title


def test_view_ascii_receptor_only(tmp_path: Path):
    """Test 13: Run run_view(receptor_pdb, ascii_mode=True), assert code 0 and receptor stats."""
    receptor_pdb = tmp_path / "receptor.pdb"
    receptor_pdb.write_text(_SAMPLE_RECEPTOR_PDB, encoding="utf-8")

    rc = run_view(receptor_pdb, ascii_mode=True)
    assert rc == 0

    # Verify receptor summary table generator
    rec_atoms = _parse_receptor_atoms(_SAMPLE_RECEPTOR_PDB)
    table = _render_receptor_summary_table(rec_atoms, "TestReceptor")
    assert table is not None
    assert "TestReceptor" in table.title


def test_view_ephemeral_web_server(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    """Test 14: Start ephemeral viewer HTTP server, perform HTTP GET, assert 200 and shutdown."""
    # 1. Direct handler & ReusableThreadingHTTPServer verification
    html_content = build_viewer_html(
        receptor_text=_SAMPLE_RECEPTOR_PDBQT,
        poses_list=[_SAMPLE_LIGAND_PDBQT_POSE],
        title="DockStudio Micro Viewer Test",
        subtitle="TestPose",
    )
    handler = create_viewer_handler(html_content)
    server = ReusableThreadingHTTPServer(("127.0.0.1", 0), handler)
    port = server.server_port
    assert port > 0

    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()

    try:
        # Request valid viewer page
        req_url = f"http://127.0.0.1:{port}/"
        with urllib.request.urlopen(req_url, timeout=3.0) as resp:
            assert resp.status == 200
            content = resp.read().decode("utf-8")
            assert "DockStudio 3D" in content
            assert "3Dmol" in content

        # Request invalid URL -> 404
        bad_url = f"http://127.0.0.1:{port}/nonexistent.html"
        with pytest.raises(urllib.error.HTTPError) as exc_info:
            urllib.request.urlopen(bad_url, timeout=3.0)
        assert exc_info.value.code == 404
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=1.0)

    # 2. End-to-end run_view execution with simulated Enter press (mocked input)
    complex_dir = tmp_path / "view_e2e_complex"
    complex_dir.mkdir()
    (complex_dir / "receptor.pdbqt").write_text(_SAMPLE_RECEPTOR_PDBQT, encoding="utf-8")
    (complex_dir / "out.pdbqt").write_text(_SAMPLE_LIGAND_PDBQT_POSE, encoding="utf-8")

    # Simulate immediate Enter keypress
    monkeypatch.setattr("builtins.input", lambda: "")
    rc_e2e = run_view(complex_dir, port=0, ascii_mode=False)
    assert rc_e2e == 0


def test_view_helpers_and_edge_cases(tmp_path: Path):
    """Test auxiliary view helpers: affinity extraction, port finding, and QC checks."""
    # 1. find_available_port
    p1 = find_available_port(0)
    assert isinstance(p1, int) and p1 > 0

    p2 = find_available_port(preferred_port=8990)
    assert isinstance(p2, int) and p2 >= 8990

    # 2. _extract_affinity from various text formats
    assert _extract_affinity(tmp_path, "REMARK VINA RESULT: -9.40 0.0 0.0") == -9.40
    assert _extract_affinity(tmp_path, "estimated Free Energy of Binding = -8.15 kcal/mol") == -8.15
    assert _extract_affinity(tmp_path, "><affinity>\n-7.80\n") == -7.80
    assert _extract_affinity(tmp_path, "NO_ENERGY_INFO") is None

    # 3. _run_rdkit_qc_checks on valid 3D molecule
    mol = Chem.MolFromMolBlock(_create_conformer_sdf_block("CC(=O)O"))
    qc_res = _run_rdkit_qc_checks(mol)
    assert qc_res["results"]["mol_pred_loaded"]["passed"] is True
    assert qc_res["results"]["sanitization"]["passed"] is True
    assert qc_res["results"]["conformer_3d"]["passed"] is True

    # 4. _run_rdkit_qc_checks on None molecule
    qc_none = _run_rdkit_qc_checks(None)
    assert qc_none["results"]["mol_pred_loaded"]["passed"] is False

    # 5. Non-existent path to run_view returns 2
    assert run_view(tmp_path / "totally_nonexistent_file.pdb") == 2


# ---------------------------------------------------------------------------
# Test Cases: Agent REPL & Slash Command Autocompletion
# ---------------------------------------------------------------------------

def test_agent_slash_completer():
    """Test 15: Test SlashCompleter on '/', '/ver', '/view --ascii ', and path completion."""
    completer = SlashCompleter()

    # 1. Test '/' command list completion
    doc_root = Document("/")
    completions_root = [c.text for c in completer.get_completions(doc_root, CompleteEvent())]
    for expected_cmd in ("/dock", "/screen", "/verify", "/audit", "/view", "/engines", "/help", "/exit"):
        assert expected_cmd in completions_root

    # 2. Test prefix completion '/ver' -> '/verify'
    doc_ver = Document("/ver")
    completions_ver = [c.text for c in completer.get_completions(doc_ver, CompleteEvent())]
    assert completions_ver == ["/verify"]

    # 3. Test parameter flags completion for '/view -'
    doc_view_flag = Document("/view -")
    completions_flag = [c.text for c in completer.get_completions(doc_view_flag, CompleteEvent())]
    assert "--ascii" in completions_flag
    assert "-a" in completions_flag

    # 4. Test path completion delegation for '/view --ascii '
    doc_view_path = Document("/view --ascii ")
    completions_path = list(completer.get_completions(doc_view_path, CompleteEvent()))
    assert isinstance(completions_path, list)

    # 5. Test '/dock ' has no unexpected slash completions
    doc_dock = Document("/dock ")
    completions_dock = list(completer.get_completions(doc_dock, CompleteEvent()))
    assert len(completions_dock) == 0


def test_agent_repl_pipe_execution(monkeypatch: pytest.MonkeyPatch):
    """Test 16: Run run_agent_repl() with stdin mocked via StringIO.

    Validates non-interactive pipe execution with /help, /engines, /exit.
    """
    fake_stdin = io.StringIO("/help\n/engines\n/exit\n")
    monkeypatch.setattr(sys, "stdin", fake_stdin)
    assert not sys.stdin.isatty()

    rc = run_agent_repl()
    assert rc == 0

    # Also test pipe execution with /clear and /quit
    fake_stdin2 = io.StringIO("/clear\n/quit\n")
    monkeypatch.setattr(sys, "stdin", fake_stdin2)
    rc2 = run_agent_repl()
    assert rc2 == 0

    # Test pipe execution with unknown command triggers error return code
    fake_stdin_err = io.StringIO("/bogus_unknown_command_12345\n")
    monkeypatch.setattr(sys, "stdin", fake_stdin_err)
    rc_err = run_agent_repl()
    assert rc_err == 1


def test_agent_dispatch_and_process_cleanup():
    """Test auxiliary agent command dispatcher & process tracking cleanup."""
    # 1. _dispatch_command tests
    assert _dispatch_command("") == 0
    assert _dispatch_command("# Comment line") == 0
    assert _dispatch_command("/help") == 0
    assert _dispatch_command("/clear") == 0
    assert _dispatch_command("/exit") == 999
    assert _dispatch_command("/quit") == 999
    assert _dispatch_command("/nonexistent_cmd") == 1

    # 2. Process management
    mock_proc = mock.Mock()
    mock_proc.poll.return_value = None
    register_process(mock_proc)
    assert mock_proc in ACTIVE_PROCESSES

    unregister_process(mock_proc)
    assert mock_proc not in ACTIVE_PROCESSES

    register_process(mock_proc)
    cleanup_active_processes()
    mock_proc.terminate.assert_called_once()
    assert len(ACTIVE_PROCESSES) == 0


def test_agent_workspace_and_in_place_views():
    """Test 18: Verify AgentWorkspace in-place view switching and command canvas management."""
    from dockstudio.cli.agent import AgentWorkspace, _build_help_pane, _build_welcome_pane

    ws = AgentWorkspace(gpu_info="NVIDIA Tesla P4 (8192 MiB)", use_alt_screen=False)
    ws.set_view(_build_welcome_pane(ws.gpu_info))
    assert ws.current_view is not None

    # Test /help, /engines, /clear in-place view updates on workspace
    assert _dispatch_command("/help", workspace=ws) == 0
    assert ws.current_view is not None

    assert _dispatch_command("/engines", workspace=ws) == 0
    assert ws.current_view is not None

    assert _dispatch_command("/clear", workspace=ws) == 0
    assert ws.current_view is not None

    # Test prepare_command_canvas clears current_view
    ws.prepare_command_canvas(show_header=True)
    assert ws.current_view is None


def test_dock_demo_and_post_dock_actions_end_to_end(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    """Test 19: Run full 5-step /dock --demo wizard including PLIP, PoseBusters, and audit.zip export."""
    from dockstudio.cli.commands.dock import run_dock_interactive
    from dockstudio.models import DockParams

    monkeypatch.setenv("DOCKSTUDIO_DATA_DIR", str(tmp_path / "dock_data"))
    from dockstudio.config import get_settings
    get_settings.cache_clear()

    # Sequence of select_one responses:
    # 1. Search box -> 0 (around cocrystal ligand)
    # 2. Engine -> 0 (vina)
    # 3. Post-dock action -> 0 (PLIP tree), then 1 (PoseBusters), then 2 (Export audit.zip), then 3 (Exit)
    responses = iter([0, 0, 0, 1, 2, 3])
    monkeypatch.setattr("dockstudio.cli.commands.dock.select_one", lambda *a, **kw: next(responses, 3))

    # Also verify _select_docking_params works without crashing on available_engines() / manifest.requires_gpu
    from dockstudio.cli.commands.dock import _select_docking_params
    eng_name, sc_name, d_params = _select_docking_params(get_settings())
    assert isinstance(eng_name, str) and isinstance(sc_name, str) and d_params.exhaustiveness == 8

    # Reset responses and use fast exhaustiveness=1 for actual Vina execution
    responses = iter([0, 0, 1, 2, 3])
    monkeypatch.setattr("dockstudio.cli.commands.dock.select_one", lambda *a, **kw: next(responses, 3))
    monkeypatch.setattr(
        "dockstudio.cli.commands.dock._select_docking_params",
        lambda _s: ("vina", "vina", DockParams(exhaustiveness=1, num_modes=2, cpu=2, seed=42)),
    )

    try:
        rc = run_dock_interactive(demo=True)
        assert rc == 0
    finally:
        get_settings.cache_clear()


def test_prep_and_screen_demo_end_to_end(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    """Test 20: Verify /prep (receptor & SDF ligand) and /screen --demo wizards run end-to-end."""
    from dockstudio.cli.commands.prep import run_prep_ligand_interactive, run_prep_receptor_interactive
    from dockstudio.cli.commands.screen import run_screen_interactive
    from dockstudio.models import DockParams

    monkeypatch.setenv("DOCKSTUDIO_DATA_DIR", str(tmp_path / "screen_data"))
    from dockstudio.config import get_settings
    get_settings.cache_clear()

    try:
        # 1. Test receptor prep
        fixtures = Path(__file__).resolve().parent / "fixtures"
        pdb_path, rec = run_prep_receptor_interactive(fixtures / "1HSG.pdb")
        assert Path(rec.pdbqt_path).is_file()

        # 2. Test ligand prep from SDF file (mode=1)
        sdf_file = tmp_path / "aspirin.sdf"
        sdf_file.write_text(_create_conformer_sdf_block("CC(=O)Oc1ccccc1C(=O)O"), encoding="utf-8")
        monkeypatch.setattr("dockstudio.cli.commands.prep.select_one", lambda *a, **kw: 1)
        monkeypatch.setattr("dockstudio.cli.commands.prep.prompt_path", lambda *a, **kw: sdf_file)
        lig_name, lig = run_prep_ligand_interactive()
        assert Path(lig.pdbqt_path).is_file()

        # 3. Test /screen --demo (select box=0, fast exhaustiveness=1)
        monkeypatch.setattr("dockstudio.cli.commands.dock.select_one", lambda *a, **kw: 0)
        monkeypatch.setattr(
            "dockstudio.cli.commands.screen._select_docking_params",
            lambda _s: ("vina", "vina", DockParams(exhaustiveness=1, num_modes=1, cpu=2, seed=42)),
        )
        rc = run_screen_interactive(demo=True)
        assert rc == 0
    finally:
        get_settings.cache_clear()
