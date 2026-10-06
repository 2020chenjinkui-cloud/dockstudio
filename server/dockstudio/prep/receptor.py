"""Receptor preparation: PDB parsing, cleaning, protonation, PDBQT via meeko."""

from __future__ import annotations

import subprocess
import sys
import uuid
from collections import defaultdict
from collections.abc import Iterable
from pathlib import Path

from openbabel import openbabel as ob

from ..logging import get_logger
from ..models import CocrystalLigand, PreparationError, PreparedReceptor

log = get_logger("prep.receptor")

_WATERS = {"HOH", "WAT", "H2O", "DOD"}
_EXCLUDE = _WATERS | {
    "NA",
    "CL",
    "K",
    "MG",
    "CA",
    "ZN",
    "MN",
    "FE",
    "CU",
    "NI",
    "CO",
    "CD",
    "HG",
    "SO4",
    "PO4",
    "GOL",
    "EDO",
    "ACT",
    "DMS",
    "MPD",
    "PEG",
    "PG4",
    "PGE",
    "FMT",
    "ACE",
    "NH4",
    "NO3",
    "BR",
    "IOD",
    "CIT",
    "TRS",
    "EPE",
    "MES",
    "BME",
    "DTT",
}


class _PdbRecord:
    def __init__(self, line: str):
        self.line = line
        self.record = line[0:6].strip()
        self.resname = line[17:20].strip()
        self.chain = line[21:22].strip() or "_"
        self.resseq = line[22:26].strip()
        self.icode = line[26:27].strip()
        self.element = line[76:78].strip().upper() if len(line) >= 78 else ""
        if not self.element:
            self.element = line[12:16].strip().lstrip("0123456789")[:1].upper()
        try:
            self.x = float(line[30:38])
            self.y = float(line[38:46])
            self.z = float(line[46:54])
        except ValueError as exc:
            raise PreparationError(f"malformed ATOM/HETATM record: {line!r}") from exc

    @property
    def key(self) -> tuple[str, str, str, str]:
        return self.resname, self.chain, self.resseq, self.icode

    @property
    def is_hydrogen(self) -> bool:
        return self.element == "H"


def _parse_pdb(text: str) -> list[_PdbRecord]:
    records: list[_PdbRecord] = []
    for line in text.splitlines():
        if len(line) >= 54 and line.startswith(("ATOM  ", "HETATM")):
            records.append(_PdbRecord(line))
    return records


def _short_id() -> str:
    return uuid.uuid4().hex[:8]


def _add_hydrogens(pdb_path: Path, out_path: Path) -> None:
    conv = ob.OBConversion()
    conv.SetInAndOutFormats("pdb", "pdb")
    mol = ob.OBMol()
    text = pdb_path.read_text(encoding="utf-8", errors="replace")
    if not conv.ReadString(mol, text):
        raise PreparationError(f"OpenBabel could not parse {pdb_path.name}")
    mol.AddHydrogens(True, True, 7.4)
    out_path.write_text(_regroup_residues(conv.WriteString(mol)), encoding="utf-8", newline="\n")


def _regroup_residues(pdb_text: str) -> str:
    """Make each residue's records contiguous (OpenBabel appends added H at the
    end; meeko requires contiguous residue blocks)."""
    atoms: list[str] = []
    other: list[str] = []
    for line in pdb_text.splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            atoms.append(line)
        else:
            other.append(line)
    order: dict[tuple[str, str, str, str], int] = {}
    for line in atoms:
        key = (line[21:22], line[22:26], line[26:27], line[17:20])
        if key not in order:
            order[key] = len(order)
    atoms.sort(key=lambda line: order[(line[21:22], line[22:26], line[26:27], line[17:20])])
    body = [line for line in other if not line.startswith(("COMPND", "AUTHOR", "END"))]
    return "\n".join(body + atoms + ["END", ""])


def _run_meeko_receptor(addh_pdb: Path, out_pdbqt: Path) -> None:
    exe = Path(sys.executable).parent / "Scripts" / "mk_prepare_receptor.exe"
    cmd: list[str] = (
        [str(exe)] if exe.is_file() else [sys.executable, "-m", "meeko.cli.mk_prepare_receptor"]
    )
    cmd += ["--read_pdb", str(addh_pdb), "-o", str(out_pdbqt.with_suffix("")), "-p", str(out_pdbqt)]
    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=600,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise PreparationError(f"meeko receptor preparation failed to run: {exc}") from exc
    if proc.returncode != 0 or not out_pdbqt.is_file():
        tail = "\n".join(((proc.stdout or "") + (proc.stderr or "")).splitlines()[-10:])
        raise PreparationError(f"meeko receptor preparation failed: {tail}")


def _obabel_pdb_to_pdbqt(pdb_path: Path, out_pdbqt: Path) -> None:
    conv = ob.OBConversion()
    conv.SetInAndOutFormats("pdb", "pdbqt")
    mol = ob.OBMol()
    text = pdb_path.read_text(encoding="utf-8", errors="replace")
    if not conv.ReadString(mol, text):
        raise PreparationError("OpenBabel could not parse PDB for fallback PDBQT conversion")
    out_pdbqt.write_text(conv.WriteString(mol), encoding="utf-8", newline="\n")


def prepare_receptor(content: bytes, filename: str, out_dir: Path) -> PreparedReceptor:
    """Prepare an uploaded PDB receptor: clean, protonate, write PDBQT."""
    if filename.lower().endswith((".cif", ".mmcif")):
        raise PreparationError("mmCIF input is not supported yet; please convert to PDB")
    if not filename.lower().endswith(".pdb"):
        raise PreparationError(f"unsupported receptor file type: {Path(filename).suffix}")

    warnings: list[str] = []
    receptor_id = _short_id()
    rec_dir = out_dir / receptor_id
    rec_dir.mkdir(parents=True, exist_ok=True)

    (rec_dir / "original.pdb").write_bytes(content)
    try:
        text = content.decode("utf-8", errors="replace")
    except Exception as exc:
        raise PreparationError(f"could not decode receptor file: {exc}") from exc
    records = _parse_pdb(text)
    if not any(r.record == "ATOM" for r in records):
        raise PreparationError("no ATOM records found in receptor PDB")

    chains = sorted({r.chain for r in records if r.record == "ATOM"})

    het_groups: dict[tuple[str, str, str, str], list[_PdbRecord]] = defaultdict(list)
    for r in records:
        if r.record == "HETATM":
            het_groups[r.key].append(r)

    cocrystals: list[CocrystalLigand] = []
    for (resname, chain, _, _), group in sorted(het_groups.items()):
        if resname in _EXCLUDE:
            continue
        heavy = [r for r in group if not r.is_hydrogen]
        if len(heavy) < 5:
            continue
        n = len(heavy)
        cocrystals.append(
            CocrystalLigand(
                resname=resname,
                chain=chain,
                n_atoms=n,
                center=(
                    round(sum(r.x for r in heavy) / n, 3),
                    round(sum(r.y for r in heavy) / n, 3),
                    round(sum(r.z for r in heavy) / n, 3),
                ),
            )
        )

    clean_lines: Iterable[str] = (
        r.line for r in records if r.record == "ATOM" or r.resname not in _WATERS
    )
    clean_path = rec_dir / "receptor_clean.pdb"
    with clean_path.open("w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(clean_lines))
        fh.write("\nEND\n")

    addh_path = rec_dir / "receptor_addH.pdb"
    pdbqt_path = rec_dir / "receptor.pdbqt"
    try:
        _add_hydrogens(clean_path, addh_path)
    except PreparationError as exc:
        warnings.append(f"protonation skipped: {exc}")
    # Docking needs a free binding pocket, so the rigid receptor PDBQT is built
    # from protein (ATOM) records only. meeko assigns charges/protonation from
    # its own residue templates; feeding it the no-hydrogen file avoids the
    # ambiguous histidine states explicit OpenBabel hydrogens can trigger.
    protein_path = rec_dir / "receptor_protein.pdb"
    with protein_path.open("w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(r.line for r in records if r.record == "ATOM"))
        fh.write("\nEND\n")
    try:
        _run_meeko_receptor(protein_path, pdbqt_path)
    except PreparationError as exc:
        log.warning("meeko receptor route failed, falling back to OpenBabel: %s", exc)
        warnings.append(f"meeko receptor preparation failed, used OpenBabel fallback: {exc}")
        _obabel_pdb_to_pdbqt(protein_path, pdbqt_path)

    return PreparedReceptor(
        id=receptor_id,
        pdbqt_path=str(pdbqt_path),
        original_path=str(rec_dir / "original.pdb"),
        chains=chains,
        cocrystal_ligands=cocrystals,
        warnings=warnings,
    )
