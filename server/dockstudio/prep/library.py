"""Ligand library preparation: bulk .sdf / .smi / .csv ingestion with TCM dual-track and strict conservation."""

from __future__ import annotations

import csv
import io
import uuid
from collections.abc import Iterator
from pathlib import Path
from typing import Any

from rdkit import Chem

from ..core.orchestrator import PrepareRegistry
from ..logging import get_logger
from ..models import PreparationError, PreparedLibrary
from .ligand import _write_ligand
from .tcm import (
    deglycosylate_tcm_mol,
    desalt_and_filter,
    generate_conformers_and_min_energy,
    scan_pains_and_tcm_alerts,
)

log = get_logger("prep.library")

MAX_LIBRARY_SIZE = 50_000


def _iter_sdf(content: str) -> Iterator[tuple[int, str, str | None]]:
    """Yield (index, mol_block, name) splitting a multi-molecule SDF on $$$$."""
    blocks = content.split("$$$$")
    for i, block in enumerate(blocks):
        block = block.strip()
        if not block:
            continue
        name = block.splitlines()[0].strip() or None
        yield i, block + "\n$$$$\n", name


def _iter_smi(content: str) -> Iterator[tuple[int, str, str | None]]:
    """Yield (index, smiles, name) from .smi lines: 'SMILES [name]'."""
    for i, line in enumerate(content.splitlines()):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split(None, 1)
        yield i, parts[0], (parts[1].strip() or None if len(parts) > 1 else None)


def _iter_csv(content: str) -> Iterator[tuple[int, str, str | None]]:
    """Yield (index, smiles, name) from CSV with a smiles column (case-insens.)."""
    reader = csv.DictReader(io.StringIO(content))
    if not reader.fieldnames:
        raise PreparationError("CSV 文件缺少表头。")
    smiles_col = next((c for c in reader.fieldnames if c.strip().lower() == "smiles"), None)
    if smiles_col is None:
        raise PreparationError(f"CSV 缺少 smiles 列（现有列: {', '.join(reader.fieldnames)}）。")
    name_col = next(
        (c for c in reader.fieldnames if c.strip().lower() in ("name", "title", "id")), None
    )
    for i, row in enumerate(reader):
        smiles = (row.get(smiles_col) or "").strip()
        if not smiles:
            continue
        yield i, smiles, ((row.get(name_col) or "").strip() or None if name_col else None)


def prepare_library(
    content: bytes,
    filename: str,
    registry: PrepareRegistry,
    name: str | None = None,
    tcm_mode: bool = False,
    min_heavy: int = 3,
    max_heavy: int = 120,
) -> PreparedLibrary:
    """Parse a .sdf/.smi/.csv upload and prepare every entry.

    Maintains strict mass conservation:
        N_raw - (N_desalt_fail + N_dedup + N_chiral_fail + N_prep_fail) = N_ready

    In TCM dual-track mode (tcm_mode=True):
    - Canonical InChIKey deduplication with source_herbs tracking.
    - Cleavable glycosides emit an aglycone derivative (<id>_aglycone) alongside
      the original glycoside (<id>_glycoside).
    """
    try:
        text = content.decode("utf-8", errors="replace")
    except Exception as exc:
        raise PreparationError(f"无法解码库文件: {exc}") from exc

    suffix = Path(filename).suffix.lower()
    if suffix == ".sdf":
        entries: list[tuple[int, str, str | None]] = list(_iter_sdf(text))
    elif suffix == ".smi":
        entries = list(_iter_smi(text))
    elif suffix == ".csv":
        entries = list(_iter_csv(text))
    else:
        raise PreparationError(f"不支持的库格式：{suffix or '<无扩展名>'}（支持 .sdf/.smi/.csv）。")

    if not entries:
        raise PreparationError("库文件中没有可解析的配体条目。")
    if len(entries) > MAX_LIBRARY_SIZE:
        raise PreparationError(
            f"库条目数 {len(entries)} 超过上限 {MAX_LIBRARY_SIZE}，请先拆分库文件。"
        )

    n_raw = len(entries)
    n_desalt_fail = 0
    n_dedup = 0
    n_chiral_fail = 0
    n_prep_fail = 0
    n_ready = 0

    library_id = uuid.uuid4().hex[:8]
    lib_dir = registry.ligands_dir / f"lib_{library_id}"
    lib_dir.mkdir(parents=True, exist_ok=True)

    ligand_ids: list[str] = []
    failed: list[dict[str, Any]] = []
    seen_inchikeys: dict[str, str] = {}
    source_herbs: dict[str, list[str]] = {}

    for index, payload, entry_name in entries:
        # A. Mol Parsing
        mol: Chem.Mol | None = None
        if suffix == ".sdf":
            mol = Chem.MolFromMolBlock(payload, sanitize=False, removeHs=False)
            if mol is not None:
                try:
                    mol.UpdatePropertyCache(strict=False)
                    Chem.FastFindRings(mol)
                    Chem.SanitizeMol(mol)
                except Exception:
                    mol = None
        else:
            payload_str = payload.strip()
            if any(ch.isspace() for ch in payload_str):
                payload_str = payload_str.split()[0]
            mol = Chem.MolFromSmiles(payload_str)

        if mol is None:
            n_desalt_fail += 1
            failed.append({
                "index": index,
                "smiles": payload[:120],
                "reason": "RDKit 无法解析分子结构",
            })
            log.warning("library entry %d failed: unparseable structure", index)
            continue

        # B. Intelligent Desalting & Whitelist/Size Filter
        clean_mol, desalt_reason = desalt_and_filter(mol, min_heavy=min_heavy, max_heavy=max_heavy)
        if clean_mol is None:
            n_desalt_fail += 1
            failed.append({
                "index": index,
                "smiles": payload[:120],
                "reason": desalt_reason or "脱盐或元素/尺寸过滤失败",
            })
            log.warning("library entry %d failed desalt/filter: %s", index, desalt_reason)
            continue

        # C. Canonical InChIKey Deduplication (in TCM mode)
        try:
            inchikey = Chem.MolToInchiKey(clean_mol)
        except Exception:
            inchikey = Chem.MolToSmiles(clean_mol)

        herb_label = entry_name or f"entry_{index}"
        if tcm_mode:
            if inchikey in seen_inchikeys:
                n_dedup += 1
                source_herbs.setdefault(inchikey, []).append(herb_label)
                log.info("Duplicate InChIKey %s merged (from %s)", inchikey, herb_label)
                continue
            seen_inchikeys[inchikey] = f"entry_{index}"
            source_herbs.setdefault(inchikey, []).append(herb_label)

        # D. Base ID & Deglycosylation Check
        base_id = uuid.uuid4().hex[:8]
        is_glycoside = False
        aglycone_mol: Chem.Mol | None = None
        if tcm_mode:
            try:
                aglycone_mol, _ = deglycosylate_tcm_mol(clean_mol)
                if aglycone_mol is not None:
                    is_glycoside = True
            except Exception as exc:
                log.warning("deglycosylation error on entry %d: %s", index, exc)

        primary_id = f"{base_id}_glycoside" if (tcm_mode and is_glycoside) else base_id
        primary_name = entry_name or primary_id

        # E. Primary Ligand Preparation
        try:
            mol_3d, min_energy = generate_conformers_and_min_energy(clean_mol)
            qc_profile = scan_pains_and_tcm_alerts(clean_mol)
            qc_profile.derivative_type = "original"
            qc_profile.is_glycoside = is_glycoside

            ligand = _write_ligand(
                mol_3d,
                source="file" if suffix == ".sdf" else "smiles",
                name=primary_name,
                out_dir=registry.ligands_dir,
                qc_profile=qc_profile,
                min_energy_kcal=min_energy,
                derivative_type="original",
                ligand_id=primary_id,
            )
            ligand_ids.append(ligand.id)
            registry.add_ligand(ligand)
            n_ready += 1
        except Exception as exc:
            n_prep_fail += 1
            failed.append({
                "index": index,
                "smiles": payload[:120],
                "reason": f"3D嵌入或PDBQT生成失败: {exc}",
            })
            log.warning("library entry %d prep failed: %s", index, exc)
            continue

        # F. Dual-Track Derived Aglycone Preparation
        if tcm_mode and is_glycoside and aglycone_mol is not None:
            aglycone_id = f"{base_id}_aglycone"
            aglycone_name = f"{primary_name}_aglycone" if entry_name else aglycone_id
            try:
                agl_3d, agl_min_energy = generate_conformers_and_min_energy(aglycone_mol)
                agl_qc = scan_pains_and_tcm_alerts(aglycone_mol)
                agl_qc.derivative_type = "aglycone"
                agl_qc.parent_glycoside_id = primary_id

                agl_ligand = _write_ligand(
                    agl_3d,
                    source="smiles",
                    name=aglycone_name,
                    out_dir=registry.ligands_dir,
                    qc_profile=agl_qc,
                    min_energy_kcal=agl_min_energy,
                    derivative_type="aglycone",
                    parent_id=primary_id,
                    ligand_id=aglycone_id,
                )
                ligand_ids.append(agl_ligand.id)
                registry.add_ligand(agl_ligand)
            except Exception as exc:
                log.warning("Failed to prepare derived aglycone for entry %d: %s", index, exc)

    # Invariant assertion
    total_losses = n_desalt_fail + n_dedup + n_chiral_fail + n_prep_fail
    assert n_raw - total_losses == n_ready, f"Loss mismatch: {n_raw} - {total_losses} != {n_ready}"

    manifest = {
        "raw": n_raw,
        "desalt": n_desalt_fail,
        "dedup": n_dedup,
        "chiral_failed": n_chiral_fail,
        "prep_failed": n_prep_fail,
        "dock_failed": 0,
        "docked_success": n_ready,
        "ready": n_ready,
        "n_ligands": len(ligand_ids),
        "total_losses": total_losses,
    }

    lib = PreparedLibrary(
        id=library_id,
        name=name or Path(filename).stem,
        path=str(lib_dir),
        ligand_ids=ligand_ids,
        failed=failed,
        n_raw=n_raw,
        n_desalt_fail=n_desalt_fail,
        n_dedup=n_dedup,
        n_chiral_fail=n_chiral_fail,
        n_prep_fail=n_prep_fail,
        n_ready=n_ready,
        source_herbs=source_herbs,
        manifest=manifest,
    )

    # Persist library.json and funnel.csv
    (lib_dir / "library.json").write_text(lib.model_dump_json(indent=2), encoding="utf-8")
    funnel_csv = (
        f"raw,{n_raw}\n"
        f"desalt,{n_desalt_fail}\n"
        f"dedup,{n_dedup}\n"
        f"chiral_failed,{n_chiral_fail}\n"
        f"prep_failed,{n_prep_fail}\n"
        f"dock_failed,0\n"
        f"docked_success,{n_ready}\n"
    )
    (lib_dir / "funnel.csv").write_text(funnel_csv, encoding="utf-8")

    log.info(
        "library %s: raw=%d, ready=%d, desalt_fail=%d, dedup=%d, prep_fail=%d, total_ligands=%d",
        library_id,
        n_raw,
        n_ready,
        n_desalt_fail,
        n_dedup,
        n_prep_fail,
        len(ligand_ids),
    )
    return lib
