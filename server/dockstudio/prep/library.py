"""Ligand library preparation: bulk .sdf / .smi / .csv ingestion with TCM dual-track and strict conservation."""

from __future__ import annotations

import csv
import io
import json
import uuid
from collections.abc import Iterator
from pathlib import Path
from typing import Any

from rdkit import Chem

from ..core.orchestrator import PrepareRegistry
from ..logging import get_logger
from ..models import PreparationError, PreparedLigand, PreparedLibrary
from .ligand import _write_ligand
from .tcm import (
    deglycosylate_tcm_mol,
    desalt_and_filter,
    enumerate_stereoisomers_tcm,
    generate_conformers_and_min_energy,
    scan_pains_and_tcm_alerts,
    standardize_tautomer_and_charge,
)

log = get_logger("prep.library")

MAX_LIBRARY_SIZE = 50_000

#: Chiral QC gate (Task Fix 1): more than this many undefined stereocentres means the
#: molecule cannot be represented faithfully and is counted as ``n_chiral_fail``.
MAX_UNSPECIFIED_STEREO_CENTERS = 3
MAX_STEREO_ISOMERS = 4


def _unassigned_centers_confined_to_aliphatic_ring(mol: Chem.Mol) -> bool:
    """True when every undefined stereocentre sits on an isolated aliphatic monocycle.

    Natural sugars / cyclitols (and the oligosaccharide chains built from them) carry a
    well-defined absolute configuration that a 2D SMILES frequently omits; enumerating
    2**n ring stereoisomers is meaningless there, so such an entry keeps a representative
    isomer instead of being rejected by the chiral QC gate.

    Every undefined centre must be a ring atom of at least one *isolated* aliphatic
    monocycle.  Acyclic polyols (``CC(O)C(O)C(O)C(O)C``), fused/bridged polycycles and
    mixed chain/ring ambiguities are still rejected as ``n_chiral_fail``.
    """
    try:
        centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    except Exception:  # pragma: no cover - defensive
        return False
    unassigned = [idx for idx, tag in centers if tag == "?"]
    if not unassigned:
        return True
    try:
        Chem.FastFindRings(mol)
        ring_info = mol.GetRingInfo()
        atom_rings = [tuple(r) for r in ring_info.AtomRings()]
    except Exception:  # pragma: no cover - defensive
        return False

    ring_membership: dict[int, int] = {}
    for ring in atom_rings:
        for atom_idx in ring:
            ring_membership[atom_idx] = ring_membership.get(atom_idx, 0) + 1

    def _is_isolated_aliphatic_ring(ring: tuple[int, ...]) -> bool:
        if any(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
            return False
        for i in range(len(ring)):
            bond = mol.GetBondBetweenAtoms(ring[i], ring[(i + 1) % len(ring)])
            if bond is None or ring_info.NumBondRings(bond.GetIdx()) > 1:
                return False  # fused / bridged
        if any(ring_membership.get(i, 0) > 1 for i in ring):
            return False  # spiro / bridgehead
        return True

    for atom_idx in unassigned:
        host_rings = [r for r in atom_rings if atom_idx in r]
        if not host_rings:
            return False  # acyclic stereocentre -> real ambiguity
        if not any(_is_isolated_aliphatic_ring(r) for r in host_rings):
            return False
    return True


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
    enumerate_stereo: bool = False,
) -> PreparedLibrary:
    """Parse a .sdf/.smi/.csv upload and prepare every entry.

    Maintains strict mass conservation:
        N_raw - (N_desalt_fail + N_dedup + N_chiral_fail + N_prep_fail) = N_ready

    In TCM dual-track mode (tcm_mode=True):
    - Tautomer/charge standardisation runs **before** InChIKey deduplication, so keto/
      enol and charge-variant forms collapse onto one canonical entry.
    - Canonical InChIKey deduplication with source_herbs tracking.
    - Cleavable glycosides emit an aglycone derivative (<id>_aglycone) alongside
      the original glycoside (<id>_glycoside).
    - Molecules with more than ``MAX_UNSPECIFIED_STEREO_CENTERS`` undefined
      stereocentres are rejected into ``n_chiral_fail`` (unless every undefined centre
      belongs to one simple aliphatic ring, e.g. a natural sugar ring).

    ``enumerate_stereo`` is opt-in: when True the additional enumerated stereoisomers
    are appended as ``derivative_type="stereoisomer"`` derivatives.  ``n_ready`` never
    counts derivatives - only the primary entry of each input row.
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
    prepared_ligands: list[PreparedLigand] = []
    failed: list[dict[str, Any]] = []
    seen_inchikeys: dict[str, str] = {}
    source_herbs: dict[str, list[str]] = {}
    n_stereo_generated = 0

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

        # B2. Tautomer / charge standardisation BEFORE InChIKey dedup (P0-4)
        try:
            clean_mol = standardize_tautomer_and_charge(clean_mol)
        except Exception as exc:
            log.warning("standardisation failed on entry %d: %s", index, exc)

        # C. Canonical InChIKey Deduplication (in TCM mode)
        try:
            inchikey = Chem.MolToInchiKey(clean_mol)
        except Exception:
            inchikey = None
        if not inchikey:
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

        # C2. Chiral QC gate: too many undefined stereocentres -> n_chiral_fail
        isomers, chiral_ok = enumerate_stereoisomers_tcm(
            clean_mol,
            max_unspecified_centers=MAX_UNSPECIFIED_STEREO_CENTERS,
            max_isomers=MAX_STEREO_ISOMERS,
        )
        if not chiral_ok and not _unassigned_centers_confined_to_aliphatic_ring(clean_mol):
            n_chiral_fail += 1
            failed.append({
                "index": index,
                "smiles": payload[:120],
                "reason": (
                    f"未定义手性中心超过 {MAX_UNSPECIFIED_STEREO_CENTERS} 个，"
                    "立体化学歧义过大已拦截"
                ),
            })
            log.warning(
                "library entry %d rejected: too many undefined stereocentres", index
            )
            continue
        if isomers:
            clean_mol = isomers[0]
        # Extra stereoisomers are opt-in derivatives and never count towards n_ready.
        extra_isomers = isomers[1:] if (enumerate_stereo and chiral_ok) else []

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
            prepared_ligands.append(ligand)
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
                prepared_ligands.append(agl_ligand)
                registry.add_ligand(agl_ligand)
            except Exception as exc:
                log.warning("Failed to prepare derived aglycone for entry %d: %s", index, exc)

        # G. Optional extra stereoisomer derivatives (opt-in, never counted in n_ready)
        for stereo_idx, stereo_mol in enumerate(extra_isomers, start=1):
            stereo_id = f"{base_id}_stereo{stereo_idx}"
            stereo_name = f"{primary_name}_stereo{stereo_idx}"
            try:
                stereo_3d, stereo_energy = generate_conformers_and_min_energy(stereo_mol)
                stereo_qc = scan_pains_and_tcm_alerts(stereo_mol)
                stereo_qc.derivative_type = "stereoisomer"
                stereo_qc.parent_glycoside_id = primary_id

                stereo_ligand = _write_ligand(
                    stereo_3d,
                    source="smiles",
                    name=stereo_name,
                    out_dir=registry.ligands_dir,
                    qc_profile=stereo_qc,
                    min_energy_kcal=stereo_energy,
                    derivative_type="stereoisomer",
                    parent_id=primary_id,
                    ligand_id=stereo_id,
                )
                ligand_ids.append(stereo_ligand.id)
                prepared_ligands.append(stereo_ligand)
                registry.add_ligand(stereo_ligand)
                n_stereo_generated += 1
            except Exception as exc:
                log.warning(
                    "Failed to prepare stereoisomer derivative for entry %d: %s", index, exc
                )

    # Invariant assertion
    total_losses = n_desalt_fail + n_dedup + n_chiral_fail + n_prep_fail
    assert n_raw - total_losses == n_ready, f"Loss mismatch: {n_raw} - {total_losses} != {n_ready}"

    # Extra dual-track derivatives (aglycones + optional stereoisomers) that entered the
    # docking pool on top of the n_ready primaries.  The screening funnel closes with
    # `n_ready + n_aglycone_generated - losses == docked_success`.
    n_aglycone_generated = len(ligand_ids) - n_ready

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
        "n_aglycone_generated": n_aglycone_generated,
        "n_stereoisomer_generated": n_stereo_generated,
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

    # Persist library.json, funnel.csv and the ligand manifest (Task Fix 5 §1)
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

    # library_manifest.json - top-level dict contract shared with screen_runner:
    #   {"ligands": [<PreparedLigand JSON>, ...], "n_raw": .., "n_ready": ..,
    #    "n_aglycone_generated": ..}
    # so ScreeningRunner can restore smiles / qc_profile / min_energy_kcal even when it
    # only receives the PreparedLibrary object.
    library_manifest = {
        "library_id": library_id,
        "name": lib.name,
        "n_raw": n_raw,
        "n_ready": n_ready,
        "n_aglycone_generated": n_aglycone_generated,
        "n_stereoisomer_generated": n_stereo_generated,
        "n_desalt_fail": n_desalt_fail,
        "n_dedup": n_dedup,
        "n_chiral_fail": n_chiral_fail,
        "n_prep_fail": n_prep_fail,
        "ligand_ids": list(ligand_ids),
        "ligands": [ligand.model_dump() for ligand in prepared_ligands],
    }
    (lib_dir / "library_manifest.json").write_text(
        json.dumps(library_manifest, ensure_ascii=False, indent=2, default=str),
        encoding="utf-8",
    )

    log.info(
        "library %s: raw=%d, ready=%d, desalt_fail=%d, dedup=%d, chiral_fail=%d, prep_fail=%d, "
        "total_ligands=%d (+%d derivatives)",
        library_id,
        n_raw,
        n_ready,
        n_desalt_fail,
        n_dedup,
        n_chiral_fail,
        n_prep_fail,
        len(ligand_ids),
        n_aglycone_generated,
    )
    return lib
