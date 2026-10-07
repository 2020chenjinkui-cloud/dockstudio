"""Receptor preparation: PDB parsing, cleaning, protonation, PDBQT via meeko."""

from __future__ import annotations

import re
import shutil
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


def _run_meeko_receptor(
    addh_pdb: Path, out_pdbqt: Path, flex_spec: str | None = None
) -> Path | None:
    """Run mk_prepare_receptor. Returns the receptor JSON path when written."""
    exe = Path(sys.executable).parent / "Scripts" / "mk_prepare_receptor.exe"
    cmd: list[str] = (
        [str(exe)] if exe.is_file() else [sys.executable, "-m", "meeko.cli.mk_prepare_receptor"]
    )
    base = out_pdbqt.with_suffix("")
    cmd += ["--read_pdb", str(addh_pdb), "-o", str(base), "-p", "--forgive_extra_bonds"]
    if flex_spec:
        cmd += ["-f", flex_spec, "-j"]
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
    # with flex, meeko emits <base>_rigid.pdbqt/_flex.pdbqt/<base>.json; without
    # flex it writes <base>.pdbqt.
    expected = base.parent / f"{base.name}_rigid.pdbqt" if flex_spec else out_pdbqt
    if proc.returncode != 0 or not expected.is_file():
        tail = "\n".join(((proc.stdout or "") + (proc.stderr or "")).splitlines()[-10:])
        raise PreparationError(f"meeko receptor preparation failed: {tail}")
    json_path = base.with_suffix(".json")
    return json_path if json_path.is_file() else None


_PDBQT_BRANCH_RE = re.compile(r"^(ROOT|ENDROOT|BRANCH|ENDBRANCH|TORSDOF|REMARK)\s*")


def _obabel_pdb_to_pdbqt(pdb_path: Path, out_pdbqt: Path) -> None:
    """OpenBabel fallback PDB->PDBQT. OB emits ROOT/BRANCH torsion tags (often
    fixed-width, e.g. "BRANCH10011002") that vina rejects in a *rigid* receptor,
    so strip them out."""
    conv = ob.OBConversion()
    conv.SetInAndOutFormats("pdb", "pdbqt")
    mol = ob.OBMol()
    text = pdb_path.read_text(encoding="utf-8", errors="replace")
    if not conv.ReadString(mol, text):
        raise PreparationError("OpenBabel could not parse PDB for fallback PDBQT conversion")
    raw = conv.WriteString(mol)
    lines = [line for line in raw.splitlines() if line.strip() and not _PDBQT_BRANCH_RE.match(line)]
    if not any(line.startswith(("ATOM", "HETATM")) for line in lines):
        raise PreparationError("OpenBabel fallback produced an empty receptor PDBQT")
    out_pdbqt.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")


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


_RES_SPEC_RE = re.compile(r"^([A-Za-z0-9_]):(\d+[A-Za-z]?)$")


def prepare_flex_variant(
    base: PreparedReceptor,
    flex_residues: list[str],
    out_dir: Path,
) -> PreparedReceptor:
    """Build a flexible-residue variant of a prepared base receptor.

    flex_residues entries are "CHAIN:RESNUM" (e.g. "A:42"), max 10.
    """
    if not flex_residues:
        raise PreparationError("flex_residues 不能为空，请至少指定一个残基。")
    if len(flex_residues) > 10:
        raise PreparationError(f"柔性残基最多 10 个（提交了 {len(flex_residues)} 个）。")

    base_dir = Path(base.original_path).parent
    protein_path = base_dir / "receptor_protein.pdb"
    if not protein_path.is_file():
        raise PreparationError(f"基础受体缺少蛋白 PDB：{protein_path}")
    available: set[tuple[str, str]] = set()
    for line in protein_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("ATOM  ") and len(line) >= 26:
            available.add((line[21:22].strip(), line[22:27].strip()))
    specs: list[str] = []
    for spec in flex_residues:
        m = _RES_SPEC_RE.match(spec.strip())
        if not m:
            raise PreparationError(
                f"残基格式错误：{spec!r}（应为 链:序号，如 A:42）。可用链: {sorted(base.chains)}"
            )
        chain, resnum = m.group(1), m.group(2)
        if chain not in base.chains:
            raise PreparationError(f"链 {chain} 不存在于受体中。可用链: {sorted(base.chains)}")
        if (chain, resnum) not in available:
            raise PreparationError(f"残基 {spec} 不存在于受体中。可用链: {sorted(base.chains)}")
        specs.append(f"{chain}:{resnum}")

    variant_id = _short_id()
    var_dir = out_dir / variant_id
    var_dir.mkdir(parents=True, exist_ok=True)
    # keep the uploaded structure traceable: link the base original
    shutil.copy2(base.original_path, var_dir / "original.pdb")
    out_base = var_dir / "receptor"
    _run_meeko_receptor(protein_path, out_base.with_suffix(".pdbqt"), flex_spec=",".join(specs))
    rigid_pdbqt = var_dir / "receptor_rigid.pdbqt"
    flex_pdbqt = var_dir / "receptor_flex.pdbqt"
    json_path = var_dir / "receptor.json"
    if not rigid_pdbqt.is_file() or not flex_pdbqt.is_file():
        raise PreparationError("meeko 柔性变体输出不完整（缺少 rigid/flex PDBQT）。")
    return PreparedReceptor(
        id=variant_id,
        pdbqt_path=str(rigid_pdbqt),
        original_path=str(var_dir / "original.pdb"),
        chains=list(base.chains),
        cocrystal_ligands=list(base.cocrystal_ligands),
        flex_residues=specs,
        flex_pdbqt_path=str(flex_pdbqt),
        receptor_json_path=str(json_path) if json_path.is_file() else None,
        parent_id=base.id,
    )
