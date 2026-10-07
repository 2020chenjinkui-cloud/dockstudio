"""Ligand library preparation: bulk .sdf / .smi / .csv ingestion."""

from __future__ import annotations

import csv
import io
import uuid
from collections.abc import Iterator
from pathlib import Path

from ..core.orchestrator import PrepareRegistry
from ..logging import get_logger
from ..models import PreparationError, PreparedLibrary
from .ligand import prepare_ligand_from_smiles

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
    content: bytes, filename: str, registry: PrepareRegistry, name: str | None = None
) -> PreparedLibrary:
    """Parse a .sdf/.smi/.csv upload and prepare every entry; failures are
    collected in `failed` without aborting the rest."""
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

    library_id = uuid.uuid4().hex[:8]
    lib_dir = registry.ligands_dir / f"lib_{library_id}"
    lib_dir.mkdir(parents=True, exist_ok=True)

    ligand_ids: list[str] = []
    failed: list[dict] = []
    for index, payload, entry_name in entries:
        try:
            if suffix == ".sdf":
                sdf_path = lib_dir / f"_entry_{index}.sdf"
                sdf_path.write_text(payload, encoding="utf-8")
                try:
                    from .ligand import prepare_ligand_from_file

                    ligand = prepare_ligand_from_file(sdf_path, registry.ligands_dir)
                finally:
                    sdf_path.unlink(missing_ok=True)
            else:
                ligand = prepare_ligand_from_smiles(payload, entry_name, registry.ligands_dir)
            ligand_ids.append(ligand.id)
            registry.add_ligand(ligand)
        except PreparationError as exc:
            failed.append({"index": index, "smiles": payload[:120], "reason": str(exc)})
            log.warning("library entry %d failed: %s", index, exc)

    lib = PreparedLibrary(
        id=library_id,
        name=name or Path(filename).stem,
        path=str(lib_dir),
        ligand_ids=ligand_ids,
        failed=failed,
    )
    (lib_dir / "library.json").write_text(lib.model_dump_json(indent=2), encoding="utf-8")
    log.info("library %s: %d prepared, %d failed", library_id, len(ligand_ids), len(failed))
    return lib
