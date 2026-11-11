"""Uni-Dock engine adapter (Linux + NVIDIA GPU; batch capable). CLI mirrors
vina with `--ligands` and `--search_mode fast|balance|detail`; output PDBQT
carries the same `REMARK VINA RESULT` rows as vina."""

from __future__ import annotations

import os
import platform
import shutil
import subprocess
from pathlib import Path
from typing import Any

from ..config import bundled_engine_binary
from ..logging import get_logger
from ..models import Box, DockingError, DockParams, DockResult, EngineUnavailableError
from .base import EngineManifest
from .vina import VinaEngine, _key_lines, annotate_error, parse_vina_out

log = get_logger("engines.unidock")


def _normalize_batch_output(matched_file: Path, lig_run_dir: Path) -> Path:
    """Archive one Uni-Dock batch output as ``<lig_run_dir>/out.pdbqt``.

    Uni-Dock's GPU batch mode dumps everything into a single shared ``--dir``
    (``<work_dir>/dock_out/<ligand_id>_out.pdbqt``), which leaves downstream
    consumers looking for ``Path(run_dir) / "out.pdbqt"`` with nothing to read
    (P1-4). Every successful ligand therefore gets the same layout as
    ``VinaEngine.dock()`` / ``BaseEngine.dock_batch()``; the raw dump in
    ``dock_out`` is intentionally **left untouched** so legacy globbing still
    works.

    Returns the path that actually holds the poses (``<lig_run_dir>/out.pdbqt``,
    or ``matched_file`` itself if the archive could not be written).
    """
    target = Path(lig_run_dir) / "out.pdbqt"
    try:
        lig_run_dir.mkdir(parents=True, exist_ok=True)
        if target.resolve() == Path(matched_file).resolve():
            return target
        if target.exists():
            target.unlink()
        try:
            os.link(matched_file, target)  # same filesystem: zero-copy
        except OSError:
            shutil.copy2(matched_file, target)
        return target
    except OSError as exc:  # pragma: no cover - permission/disk issues
        log.warning("unidock batch: could not archive %s as %s: %s", matched_file, target, exc)
        return Path(matched_file)


class UniDockEngine(VinaEngine):
    """Uni-Dock: GPU-batched molecular docking."""

    manifest = EngineManifest(
        name="unidock",
        scoring_functions=["vina", "vinardo", "ad4"],
        supports_flex_residues=False,
        supports_batch_ligands=True,
        requires_gpu=True,
        platforms=["linux"],
        homepage="https://github.com/dptech-corp/Uni-Dock",
    )

    def _binary(self) -> Path:
        if platform.system() != "Linux":
            raise EngineUnavailableError(
                "Uni-Dock 需要 Linux + NVIDIA GPU，Windows 用户请用 WSL2"
                "（指引见 docs/deploy.md#wsl2）。"
            )
        env = os.environ.get("DOCKSTUDIO_UNIDOCK")
        if env and Path(env).is_file():
            return Path(env)
        found = shutil.which("unidock")
        if found:
            return Path(found)
        bundled = bundled_engine_binary("unidock")
        if bundled is not None:
            return bundled
        raise EngineUnavailableError(
            "Uni-Dock 二进制未找到：设置 DOCKSTUDIO_UNIDOCK 或放入 resources/bin/linux-x86_64。"
        )

    def _gpu_ok(self) -> tuple[bool, str]:
        if shutil.which("nvidia-smi"):
            try:
                proc = subprocess.run(
                    ["nvidia-smi", "-L"], capture_output=True, timeout=15, check=False
                )
                if proc.returncode == 0 and b"GPU" in (proc.stdout or b""):
                    return True, ""
            except (OSError, subprocess.TimeoutExpired):
                pass
        for candidate in ("libcuda.so", "libcuda.so.1"):
            for libdir in ("/usr/lib/x86_64-linux-gnu", "/usr/lib64", "/usr/local/cuda/lib64"):
                if (Path(libdir) / candidate).exists():
                    return True, ""
        return False, "检测到 Uni-Dock 但无 NVIDIA GPU（找不到 nvidia-smi 或 libcuda）"

    def is_available(self) -> tuple[bool, str | None]:
        try:
            binary = self._binary()
            self.manifest.version = self._version(binary)
        except EngineUnavailableError as exc:
            return False, str(exc)
        ok, reason = self._gpu_ok()
        if not ok:
            return False, reason
        return True, None

    def _build_cmd(
        self,
        binary: str,
        box: Box,
        params: DockParams,
        receptor_name: str,
        ligand_name: str,
        out_name: str,
    ) -> list[str]:
        # Uni-Dock: --ligands (batch spelling), --search_mode in place of
        # --exhaustiveness, no --cpu/--energy_range.
        mode = (
            "fast"
            if params.exhaustiveness <= 4
            else ("balance" if params.exhaustiveness <= 16 else "detail")
        )
        return [
            binary,
            "--receptor",
            receptor_name,
            "--ligands",
            ligand_name,
            "--center_x",
            str(box.center[0]),
            "--center_y",
            str(box.center[1]),
            "--center_z",
            str(box.center[2]),
            "--size_x",
            str(box.size[0]),
            "--size_y",
            str(box.size[1]),
            "--size_z",
            str(box.size[2]),
            "--search_mode",
            mode,
            "--num_modes",
            str(params.num_modes),
            "--scoring",
            params.scoring,
            "--out",
            out_name,
        ]

    def _build_batch_cmd(
        self,
        binary: str,
        box: Box,
        params: DockParams,
        receptor_path: str,
        ligand_index_path: str,
        out_dir_path: str,
    ) -> list[str]:
        """Constructs the Uni-Dock CLI command for multi-ligand GPU batch docking."""
        mode = (
            "fast"
            if params.exhaustiveness <= 4
            else ("balance" if params.exhaustiveness <= 16 else "detail")
        )
        cmd = [
            binary,
            "--receptor",
            receptor_path,
            "--ligand_index",
            ligand_index_path,
            "--dir",
            out_dir_path,
            "--center_x",
            str(box.center[0]),
            "--center_y",
            str(box.center[1]),
            "--center_z",
            str(box.center[2]),
            "--size_x",
            str(box.size[0]),
            "--size_y",
            str(box.size[1]),
            "--size_z",
            str(box.size[2]),
            "--search_mode",
            mode,
            "--num_modes",
            str(params.num_modes),
            "--scoring",
            params.scoring,
        ]
        if params.seed is not None:
            cmd += ["--seed", str(params.seed)]
        return cmd

    def dock_batch(
        self,
        receptor_pdbqt: Path,
        batch: list[tuple[str, Path]] | list[Path],
        box: Box,
        params: DockParams,
        work_dir: Path,
        flex_pdbqt: Path | None = None,
        **kwargs: Any,
    ) -> dict[str, DockResult]:
        """Native GPU-accelerated batch docking using Uni-Dock --ligand_index and --dir.

        Executes the entire batch in a single GPU process invocation, avoiding repeated
        receptor grid recalculations and process spawning overhead.
        """
        # 1. Check engine availability
        binary = self._binary()
        version = self._version(binary)
        self.manifest.version = version

        # 2. Argument-order polymorphism support
        if isinstance(receptor_pdbqt, (list, tuple)) and not isinstance(batch, (list, tuple)):
            raw_batch = receptor_pdbqt
            receptor_path = Path(batch)
        else:
            raw_batch = batch
            receptor_path = Path(receptor_pdbqt)

        # 3. Normalize batch items to list of (ligand_id, Path)
        normalized: list[tuple[str, Path]] = []
        for item in raw_batch:
            if isinstance(item, tuple) and len(item) == 2:
                lig_id, p = item
                normalized.append((str(lig_id), Path(p)))
            elif isinstance(item, (str, Path)):
                p = Path(item)
                normalized.append((p.stem, p))
            elif hasattr(item, "id") and hasattr(item, "pdbqt_path"):
                normalized.append((str(item.id), Path(item.pdbqt_path)))

        if not normalized:
            return {}

        work_dir = Path(work_dir)
        work_dir.mkdir(parents=True, exist_ok=True)
        staging_dir = work_dir / "staging_ligands"
        staging_dir.mkdir(parents=True, exist_ok=True)
        out_dir = work_dir / "dock_out"
        out_dir.mkdir(parents=True, exist_ok=True)

        # 4. Stage receptor locally in work_dir
        local_receptor = work_dir / "receptor.pdbqt"
        if receptor_path.resolve() != local_receptor.resolve():
            shutil.copy2(receptor_path, local_receptor)

        # 5. Disambiguate ligand filenames into staging_dir
        # Since ligands in DockStudio are often named 'ligand.pdbqt' in separate folders,
        # staging them as '<ligand_id>.pdbqt' ensures Uni-Dock produces '<ligand_id>_out.pdbqt'
        # without filename collisions.
        valid_items: list[tuple[str, Path]] = []
        for lig_id, lig_path in normalized:
            if not lig_path.is_file():
                log.warning("UniDock batch: ligand file %s not found for %s", lig_path, lig_id)
                continue
            staged_path = staging_dir / f"{lig_id}.pdbqt"
            if not staged_path.exists():
                try:
                    os.symlink(lig_path.resolve(), staged_path)
                except OSError:
                    try:
                        os.link(lig_path.resolve(), staged_path)
                    except OSError:
                        shutil.copy2(lig_path, staged_path)
            valid_items.append((lig_id, staged_path))

        if not valid_items:
            log.warning("UniDock batch: no valid ligand files found in batch of %d", len(normalized))
            return {}

        # 6. Write --ligand_index file
        index_file = work_dir / "ligand_index.txt"
        index_lines = [str(staged_path.absolute()) for _, staged_path in valid_items]
        index_file.write_text("\n".join(index_lines) + "\n", encoding="utf-8")

        # 7. Build command
        cmd = self._build_batch_cmd(
            binary=str(binary),
            box=box,
            params=params,
            receptor_path=str(local_receptor.resolve()),
            ligand_index_path=str(index_file.resolve()),
            out_dir_path=str(out_dir.resolve()),
        )

        log.info("running unidock batch (size=%d): %s", len(valid_items), " ".join(cmd[:8]) + "...")

        # 8. Launch batch subprocess
        try:
            proc = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                encoding="utf-8",
                errors="replace",
                cwd=work_dir,
                env=self._subprocess_env(),
                creationflags=subprocess.CREATE_NEW_PROCESS_GROUP
                if platform.system() == "Windows"
                else 0,
                start_new_session=platform.system() != "Windows",
            )
        except OSError as exc:
            raise EngineUnavailableError(f"failed to launch unidock: {exc}") from exc

        self._register_proc(proc)

        # Timeout scaling for large batches
        batch_timeout = max(self.settings.job_timeout_s, len(valid_items) * 10)
        try:
            stdout, stderr = proc.communicate(timeout=batch_timeout)
        except subprocess.TimeoutExpired as exc:
            # Drop it from the table first so the same pid is never killed twice,
            # then terminate only this batch process (sibling workers keep going).
            self._unregister_proc(proc)
            self._terminate_proc(proc)
            raise DockingError(
                f"unidock batch timed out after {batch_timeout}s for {len(valid_items)} ligands"
            ) from exc
        finally:
            self._unregister_proc(proc)

        run_log = (stdout or "") + "\n--- stderr ---\n" + (stderr or "")

        # 9. Parse output files from out_dir and normalize each ligand's output
        #    to the standard `<run_dir>/out.pdbqt` layout (P1-4).
        results: dict[str, DockResult] = {}
        for lig_id, _ in valid_items:
            # Candidate output paths in order of preference
            candidate_paths = [
                out_dir / f"{lig_id}_out.pdbqt",
                out_dir / f"{lig_id}.pdbqt",
            ]
            matched_file: Path | None = None
            for cp in candidate_paths:
                if cp.is_file():
                    matched_file = cp
                    break
            if matched_file is None:
                # Fallback glob check in out_dir
                globs = list(out_dir.glob(f"{lig_id}*out*.pdbqt"))
                if globs:
                    matched_file = globs[0]

            if matched_file is None or not matched_file.is_file():
                log.warning("unidock batch: output file missing for ligand %s", lig_id)
                continue

            lig_run_dir = work_dir / f"dock_{lig_id}"
            out_pose = _normalize_batch_output(matched_file, lig_run_dir)

            try:
                poses = parse_vina_out(out_pose)
                results[lig_id] = DockResult(
                    engine=self.manifest.name,
                    engine_version=version,
                    params=params,
                    box=box,
                    poses=poses,
                    log=f"Uni-Dock GPU batch output for {lig_id}",
                    run_dir=str(lig_run_dir),
                )
            except Exception as exc:
                log.warning("unidock batch: failed to parse poses for ligand %s: %s", lig_id, exc)
                continue

        # 10. Check overall exit status
        if proc.returncode != 0 and len(results) == 0:
            hint = annotate_error(run_log)
            msg = f"unidock batch exited with code {proc.returncode}: {_key_lines(run_log)}"
            raise DockingError(f"{msg}\n{hint}" if hint else msg)
        elif proc.returncode != 0:
            log.warning(
                "unidock batch returned exit code %d, but %d/%d ligands succeeded",
                proc.returncode,
                len(results),
                len(valid_items),
            )

        return results
