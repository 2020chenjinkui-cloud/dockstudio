"""DockStudio CLI: Micro 3D Web Viewer & Terminal ASCII Fallback (/view).

Provides an ephemeral zero-dependency ThreadingHTTPServer serving an interactive
3Dmol.js viewer with dynamic port negotiation, SSH tunneling instructions, and
immediate cleanup upon Enter or Ctrl+C. Supports terminal --ascii fallback with
genuine Euclidean contact analysis, RDKit property inspection, and dynamic QC.
"""

from __future__ import annotations

import getpass
import json
import math
import re
import socket
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from rich.table import Table
from rich.text import Text

from ...analysis.export import _split_models
from ...config import get_settings
from ...core.project import ProjectStore
from ..ui import (
    COLOR_ERROR,
    console,
    make_header,
    make_pane,
    render_interaction_tree,
    render_posecheck_table,
)

# Single-Page 3Dmol.js HTML Template with Claude Code Dark Theme
_VIEWER_HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="zh-CN">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>DockStudio 3D Micro Viewer - __TITLE__</title>
  <script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.0.4/3Dmol-min.js"></script>
  <style>
    * { box-sizing: border-box; margin: 0; padding: 0; }
    body {
      background-color: #18181B;
      color: #F4F4F5;
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, monospace;
      overflow: hidden;
      width: 100vw;
      height: 100vh;
    }
    #viewport {
      width: 100%;
      height: 100%;
      position: absolute;
      top: 0;
      left: 0;
      z-index: 1;
    }
    .hud-panel {
      position: absolute;
      z-index: 10;
      background: rgba(39, 39, 42, 0.88);
      backdrop-filter: blur(12px);
      border: 1px solid #3F3F46;
      border-radius: 8px;
      padding: 12px 16px;
      box-shadow: 0 8px 32px rgba(0, 0, 0, 0.4);
    }
    .top-bar {
      top: 16px;
      left: 16px;
      display: flex;
      align-items: center;
      gap: 12px;
    }
    .brand {
      color: #D77757;
      font-weight: bold;
      font-size: 15px;
      letter-spacing: 0.5px;
    }
    .title-badge {
      background: #2C323E;
      color: #93A5FF;
      font-size: 12px;
      padding: 3px 8px;
      border-radius: 4px;
      font-weight: 500;
    }
    .controls-panel {
      bottom: 24px;
      left: 50%;
      transform: translateX(-50%);
      display: flex;
      align-items: center;
      gap: 14px;
      flex-wrap: wrap;
    }
    select, button {
      background: #27272A;
      color: #F4F4F5;
      border: 1px solid #52525B;
      padding: 6px 12px;
      border-radius: 6px;
      font-size: 13px;
      cursor: pointer;
      outline: none;
      transition: all 0.2s;
    }
    select:hover, button:hover {
      border-color: #D77757;
      background: #3F3F46;
    }
    button.active {
      background: #D77757;
      border-color: #D77757;
      color: #FFFFFF;
    }
    .status-text {
      font-size: 12px;
      color: #A1A1AA;
    }
    .key-hint {
      color: #71717A;
      font-size: 11px;
    }
  </style>
</head>
<body>
  <div id="viewport"></div>

  <div class="hud-panel top-bar">
    <span class="brand">DockStudio 3D</span>
    <span class="title-badge">__SUBTITLE__</span>
    <span class="status-text" id="status-display">正在初始化三维模型...</span>
  </div>

  <div class="hud-panel controls-panel">
    <label for="poseSelect" class="status-text">构象选择:</label>
    <select id="poseSelect" onchange="switchPose(parseInt(this.value))">
      __POSE_OPTIONS__
    </select>

    <button id="btnPocket" class="active" onclick="togglePocket()">5Å 口袋残基</button>
    <button id="btnCartoon" class="active" onclick="toggleCartoon()">受体骨架</button>
    <button id="btnSurface" onclick="toggleSurface()">结合袋表面</button>
    <button onclick="recenter()">聚焦配体</button>
    <button onclick="toggleSpin()">旋转动画</button>
  </div>

  <!-- Raw Molecular Data Embeds -->
  <script type="text/plain" id="receptor-raw">__RECEPTOR_DATA__</script>
  <script type="application/json" id="poses-raw">__POSES_JSON__</script>

  <script>
    let viewer = null;
    let receptorModel = null;
    let poseModels = [];
    let currentPoseIdx = 0;
    let showPocket = true;
    let showCartoon = true;
    let showSurface = false;
    let surfaceObj = null;
    let isSpinning = false;

    window.addEventListener('DOMContentLoaded', () => {
      const container = document.getElementById('viewport');
      viewer = $3Dmol.createViewer(container, { backgroundColor: '#18181B' });

      // 1. Load Receptor
      const recData = document.getElementById('receptor-raw').textContent.trim();
      if (recData) {
        const fmt = recData.includes('ATOM') ? 'pdb' : 'pdbqt';
        receptorModel = viewer.addModel(recData, fmt);
        viewer.setStyle({ model: receptorModel }, { cartoon: { color: 'spectrum', opacity: 0.85 } });
      }

      // 2. Load Poses
      const rawPoses = JSON.parse(document.getElementById('poses-raw').textContent || '[]');
      rawPoses.forEach((pText, idx) => {
        const m = viewer.addModel(pText, 'pdbqt');
        viewer.setStyle({ model: m }, { stick: { colorscheme: 'greenCarbon', radius: 0.22 } });
        m.setShow(idx === 0);
        poseModels.push(m);
      });

      updateView();
      recenter();
      document.getElementById('status-display').textContent = `已载入: ${poseModels.length} 个构象`;
    });

    function updateView() {
      if (!viewer) return;
      if (receptorModel) {
        if (showCartoon) {
          viewer.setStyle({ model: receptorModel }, { cartoon: { color: 'spectrum', opacity: 0.85 } });
        } else {
          viewer.setStyle({ model: receptorModel }, {});
        }

        if (showPocket && poseModels.length > 0) {
          const activePose = poseModels[currentPoseIdx];
          viewer.addStyle(
            { model: receptorModel, within: { distance: 5.0, sel: { model: activePose } } },
            { stick: { colorscheme: 'amino', radius: 0.16 } }
          );
        }
      }
      viewer.render();
    }

    function switchPose(idx) {
      currentPoseIdx = idx;
      poseModels.forEach((m, i) => m.setShow(i === idx));
      updateView();
      recenter();
    }

    function togglePocket() {
      showPocket = !showPocket;
      document.getElementById('btnPocket').classList.toggle('active', showPocket);
      updateView();
    }

    function toggleCartoon() {
      showCartoon = !showCartoon;
      document.getElementById('btnCartoon').classList.toggle('active', showCartoon);
      updateView();
    }

    function toggleSurface() {
      showSurface = !showSurface;
      document.getElementById('btnSurface').classList.toggle('active', showSurface);
      if (showSurface && receptorModel) {
        surfaceObj = viewer.addSurface($3Dmol.SurfaceType.VDW, { opacity: 0.5, color: 'white' }, { model: receptorModel });
      } else if (surfaceObj) {
        viewer.removeSurface(surfaceObj);
        surfaceObj = null;
      }
      viewer.render();
    }

    function recenter() {
      if (poseModels.length > 0 && poseModels[currentPoseIdx]) {
        viewer.zoomTo({ model: poseModels[currentPoseIdx] }, 1000);
      } else {
        viewer.zoomTo();
      }
    }

    function toggleSpin() {
      isSpinning = !isSpinning;
      viewer.spin(isSpinning);
    }
  </script>
</body>
</html>
"""


class ReusableThreadingHTTPServer(ThreadingHTTPServer):
    allow_reuse_address = True
    daemon_threads = True


def find_available_port(preferred_port: int = 8080, max_attempts: int = 25) -> int:
    """Finds an available local port starting from preferred_port or OS assignment."""
    if preferred_port <= 0:
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            s.bind(("127.0.0.1", 0))
            return s.getsockname()[1]

    for p in range(preferred_port, preferred_port + max_attempts):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
            try:
                s.bind(("127.0.0.1", p))
                return p
            except OSError:
                continue

    # Fallback to ephemeral port 0
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


def build_viewer_html(
    receptor_text: str,
    poses_list: list[str],
    title: str = "分子对接构象预览",
    subtitle: str = "3D 相互作用",
) -> str:
    """Generates the full self-contained HTML page embedding structure data."""
    options = []
    for i in range(len(poses_list)):
        options.append(f'<option value="{i}">Pose #{i + 1}</option>')
    options_html = "\n      ".join(options) if options else '<option value="0">默认构象</option>'

    html = _VIEWER_HTML_TEMPLATE
    html = html.replace("__TITLE__", title)
    html = html.replace("__SUBTITLE__", subtitle)
    html = html.replace("__POSE_OPTIONS__", options_html)
    html = html.replace("__RECEPTOR_DATA__", receptor_text)
    html = html.replace("__POSES_JSON__", json.dumps(poses_list, ensure_ascii=False))
    return html


def create_viewer_handler(html_content: str):
    class ViewerHandler(BaseHTTPRequestHandler):
        def do_GET(self):
            if self.path in ("/", "/index.html"):
                encoded = html_content.encode("utf-8")
                self.send_response(200)
                self.send_header("Content-Type", "text/html; charset=utf-8")
                self.send_header("Content-Length", str(len(encoded)))
                self.send_header("Connection", "close")
                self.end_headers()
                self.wfile.write(encoded)
            else:
                self.send_response(404)
                self.end_headers()

        def log_message(self, format, *args):
            # Suppress default noisy access logs to keep terminal pristine
            return

    return ViewerHandler


# ---------------------------------------------------------------------------
# Scientific Property Extraction & Dynamic Interaction Helpers
# ---------------------------------------------------------------------------


def _extract_affinity(target_path: Path, pose_text: str = "") -> float | None:
    """Dynamically parses docking affinity (kcal/mol) from pose text or job metadata."""
    patterns = [
        r"VINA RESULT:\s*([-\d.]+)",
        r"estimated Free Energy of Binding\s*=\s*([-\d.]+)",
        r">\s*<[^>]*(?:affinity|docking_score|score|energy)[^>]*>\s*\r?\n\s*([-\d.]+)",
        r"REMARK\s+.*(?:affinity|docking_score)[:=\s]+([-\d.]+)",
    ]
    for pat in patterns:
        m = re.search(pat, pose_text, re.IGNORECASE)
        if m:
            try:
                return float(m.group(1))
            except ValueError:
                pass

    if target_path.is_dir():
        job_file = target_path / "job.json"
        if job_file.is_file():
            try:
                data = json.loads(job_file.read_text(encoding="utf-8"))
                poses = data.get("result", {}).get("poses", [])
                if poses and "affinity" in poses[0]:
                    return float(poses[0]["affinity"])
            except Exception:
                pass
        for p_json in [target_path / "analysis" / "pose_1.json", target_path / "analysis_1.json"]:
            if p_json.is_file():
                try:
                    data = json.loads(p_json.read_text(encoding="utf-8"))
                    if "affinity" in data:
                        return float(data["affinity"])
                except Exception:
                    pass
    elif target_path.is_file():
        job_file = target_path.parent / "job.json"
        if job_file.is_file():
            try:
                data = json.loads(job_file.read_text(encoding="utf-8"))
                poses = data.get("result", {}).get("poses", [])
                if poses and "affinity" in poses[0]:
                    return float(poses[0]["affinity"])
            except Exception:
                pass
    return None


def _parse_receptor_atoms(receptor_text: str) -> list[dict[str, Any]]:
    """Extracts heavy atoms and residue labels from receptor PDB/PDBQT text."""
    atoms: list[dict[str, Any]] = []
    for line in receptor_text.splitlines():
        if line.startswith(("ATOM  ", "HETATM")) and len(line) >= 54:
            elem = line[76:78].strip().upper() if len(line) >= 78 else ""
            if not elem:
                elem = line[12:16].strip().lstrip("0123456789")[:1].upper()
            if elem == "H":
                continue
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            resname = line[17:20].strip()
            chain = line[21:22].strip() or "A"
            resnr = line[22:26].strip()
            atoms.append({
                "elem": elem,
                "pos": (x, y, z),
                "resname": resname,
                "chain": chain,
                "resnr": resnr,
                "name": line[12:16].strip(),
                "res_label": f"{resname}{resnr}:{chain}",
            })
    return atoms


def _parse_ligand_atoms_and_mol(pose_text: str) -> tuple[list[dict[str, Any]], Chem.Mol | None]:
    """Parses ligand atomic 3D coordinates and builds an RDKit Mol object."""
    mol: Chem.Mol | None = None
    atoms: list[dict[str, Any]] = []

    # 1. Try SDF MolBlock
    if "$$$$" in pose_text or "M  END" in pose_text:
        try:
            mol = Chem.MolFromMolBlock(pose_text, sanitize=True, removeHs=False)
        except Exception:
            mol = None
        if mol is None:
            try:
                mol = Chem.MolFromMolBlock(pose_text, sanitize=False, removeHs=False)
            except Exception:
                mol = None

    # 2. Try OpenBabel PDBQT -> SDF conversion to preserve 3D conformer
    if mol is None and ("ATOM" in pose_text or "HETATM" in pose_text):
        try:
            from openbabel import pybel
            ob = pybel.readstring("pdbqt", pose_text)
            mol = Chem.MolFromMolBlock(ob.write("sdf"), sanitize=False, removeHs=False)
        except Exception:
            mol = None

    # 3. Try PDB Block
    if mol is None and ("ATOM" in pose_text or "HETATM" in pose_text):
        try:
            mol = Chem.MolFromPDBBlock(pose_text, sanitize=False, removeHs=False)
        except Exception:
            mol = None

    # 4. Fallback: Parse REMARK SMILES for 2D topological descriptors
    if mol is None:
        m_smi = re.search(r"REMARK\s+SMILES\s+([^\r\n]+)", pose_text)
        if m_smi:
            try:
                mol = Chem.MolFromSmiles(m_smi.group(1).strip())
            except Exception:
                mol = None

    # Extract 3D atom positions
    if mol is not None and mol.GetNumConformers() > 0:
        conf = mol.GetConformer()
        for i in range(mol.GetNumAtoms()):
            a = mol.GetAtomWithIdx(i)
            elem = a.GetSymbol().upper()
            if elem == "H":
                continue
            pt = conf.GetAtomPosition(i)
            atoms.append({"elem": elem, "pos": (pt.x, pt.y, pt.z), "idx": i})
    else:
        # Fallback: extract heavy atoms directly from ATOM/HETATM lines
        for line in pose_text.splitlines():
            if line.startswith("ENDMDL"):
                break
            if line.startswith(("ATOM  ", "HETATM")) and len(line) >= 54:
                elem = line[76:78].strip().upper() if len(line) >= 78 else ""
                if not elem:
                    elem = line[12:16].strip().lstrip("0123456789")[:1].upper()
                if elem == "H":
                    continue
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    atoms.append({"elem": elem, "pos": (x, y, z), "idx": len(atoms)})
                except ValueError:
                    continue

    return atoms, mol


def _calculate_dynamic_interactions(
    rec_atoms: list[dict[str, Any]],
    lig_atoms: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Calculates genuine Euclidean close contacts and hydrogen bonds."""
    if not rec_atoms or not lig_atoms:
        return []

    hbonds: dict[str, dict[str, Any]] = {}
    hydrophobics: dict[str, dict[str, Any]] = {}

    polar_elems = {"N", "O"}
    hydrophobic_res = {
        "ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "TYR", "HIS", "CYS"
    }

    for la in lig_atoms:
        lp = la["pos"]
        le = la["elem"]
        for ra in rec_atoms:
            rp = ra["pos"]
            re = ra["elem"]
            d = math.dist(lp, rp)
            res_label = ra["res_label"]

            # Hydrogen bonds: polar N/O donor/acceptor <= 3.5 Å
            if d <= 3.5 and le in polar_elems and re in polar_elems:
                if res_label not in hbonds or d < hbonds[res_label]["distance"]:
                    hbonds[res_label] = {
                        "type": "hbond",
                        "residue": res_label,
                        "distance": round(d, 2),
                        "detail": f"{ra['name']}-{le} 极性氢键接触",
                    }
            # Hydrophobic / close contacts <= 4.0 Å
            elif d <= 4.0:
                is_hydro = (le == "C" and re == "C" and ra["resname"] in hydrophobic_res)
                if is_hydro:
                    if res_label not in hydrophobics or d < hydrophobics[res_label]["distance"]:
                        hydrophobics[res_label] = {
                            "type": "hydrophobic",
                            "residue": res_label,
                            "distance": round(d, 2),
                            "detail": f"{ra['resname']} 侧链疏水接触",
                        }
                elif d <= 3.8 and res_label not in hbonds and res_label not in hydrophobics:
                    hydrophobics[res_label] = {
                        "type": "hydrophobic",
                        "residue": res_label,
                        "distance": round(d, 2),
                        "detail": f"近距离空间范德华接触 ({round(d, 2)} Å)",
                    }

    items = list(hbonds.values()) + list(hydrophobics.values())
    items.sort(key=lambda x: (0 if x["type"] == "hbond" else 1, x.get("distance", 99.0)))
    return items


def _render_ligand_properties_table(mol: Chem.Mol, ligand_name: str) -> Table:
    """Renders authentic RDKit physicochemical descriptors for single ligand view."""
    table = Table(
        title=f"[bold #D77757]配体物理化学性质概况 (RDKit 计算) ── {ligand_name}[/bold #D77757]",
        border_style="#52525B",
        header_style="bold #FFFFFF",
    )
    table.add_column("属性指标 (Property)", style="bold #FFFFFF", width=26)
    table.add_column("数值 (Value)", style="#4EBA65", justify="right", width=16)
    table.add_column("参考范围 / 药效团特征 (Rule of 5 / Veber)", style="#999999")

    formula = rdMolDescriptors.CalcMolFormula(mol)
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    tpsa = Descriptors.TPSA(mol)
    rot_bonds = Descriptors.NumRotatableBonds(mol)
    heavy_atoms = mol.GetNumHeavyAtoms()
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)
    aromatic_rings = Descriptors.NumAromaticRings(mol)

    table.add_row("分子式 (Formula)", formula, "元素组成守恒")
    table.add_row("分子量 (MW)", f"{mw:.2f} Da", "[#4EBA65]符合[/#4EBA65]" if mw <= 500 else "[#FFC107]MW > 500 Da[/#FFC107]")
    table.add_row("脂水分配系数 (cLogP)", f"{logp:.2f}", "[#4EBA65]符合[/#4EBA65]" if logp <= 5.0 else "[#FFC107]LogP > 5.0[/#FFC107]")
    table.add_row("极性表面积 (TPSA)", f"{tpsa:.2f} Å²", "[#4EBA65]良好透膜性[/#4EBA65]" if tpsa <= 140 else "[#FFC107]TPSA > 140 Å²[/#FFC107]")
    table.add_row("可旋转键数 (Rotatable Bonds)", str(rot_bonds), "[#4EBA65]柔性合规[/#4EBA65]" if rot_bonds <= 10 else "[#FFC107]柔性较高 (>10)[/#FFC107]")
    table.add_row("重原子数 (Heavy Atoms)", str(heavy_atoms), "非氢骨架原子总数")
    table.add_row("氢键供体数 (HBD)", str(hbd), "[#4EBA65]符合[/#4EBA65]" if hbd <= 5 else "[#FFC107]HBD > 5[/#FFC107]")
    table.add_row("氢键受体数 (HBA)", str(hba), "[#4EBA65]符合[/#4EBA65]" if hba <= 10 else "[#FFC107]HBA > 10[/#FFC107]")
    table.add_row("芳香环数量 (Aromatic Rings)", str(aromatic_rings), "芳香体系环数")
    return table


def _render_receptor_summary_table(rec_atoms: list[dict[str, Any]], receptor_name: str) -> Table:
    """Renders structural statistics for receptor-only view."""
    table = Table(
        title=f"[bold #D77757]受体蛋白结构概况 ── {receptor_name}[/bold #D77757]",
        border_style="#52525B",
        header_style="bold #FFFFFF",
    )
    table.add_column("结构指标 (Feature)", style="bold #FFFFFF", width=26)
    table.add_column("数值 / 统计 (Value)", style="#93A5FF", justify="left", width=22)
    table.add_column("说明 (Description)", style="#999999")

    chains = sorted(list({a["chain"] for a in rec_atoms}))
    res_set = {a["res_label"] for a in rec_atoms}
    res_types = sorted(list({a["resname"] for a in rec_atoms}))
    heavy_count = len(rec_atoms)

    table.add_row("链标识 (Chains)", ", ".join(chains) if chains else "A", f"共 {len(chains)} 条肽链")
    table.add_row("残基总数 (Residues)", str(len(res_set)), "非氢标准残基与辅因子计数")
    table.add_row("重原子总数 (Heavy Atoms)", str(heavy_count), "蛋白质非氢骨架与侧链原子")
    table.add_row("主要残基种类 (Residue Types)", f"{len(res_types)} 种", ", ".join(res_types[:8]) + ("..." if len(res_types) > 8 else ""))
    return table


def _run_rdkit_qc_checks(
    mol: Chem.Mol | None,
    rec_atoms: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    """Runs genuine RDKit sanitization and quality checks replacing static mock."""
    if mol is None:
        return {
            "results": {
                "mol_pred_loaded": {"passed": False, "message": "无法解析分子拓扑结构"},
                "sanitization": {"passed": False, "message": "跳过 (分子解析失败)"},
                "conformer_3d": {"passed": False, "message": "跳过 (分子解析失败)"},
            }
        }

    results: dict[str, dict[str, Any]] = {}
    results["mol_pred_loaded"] = {"passed": True, "message": "分子拓扑正常载入"}

    # 1. Sanitization: Valence & aromaticity
    try:
        mol_copy = Chem.Mol(mol)
        Chem.SanitizeMol(mol_copy)
        results["sanitization"] = {"passed": True, "message": "化合价与芳香性有效"}
    except Exception as exc:
        results["sanitization"] = {"passed": False, "message": f"化合价/芳香性异常: {exc}"}

    # 2. 3D Conformer check
    has_3d = False
    if mol.GetNumConformers() > 0:
        conf = mol.GetConformer()
        if conf.Is3D():
            zs = [conf.GetAtomPosition(i).z for i in range(mol.GetNumAtoms())]
            if max(zs) - min(zs) > 0.05 or mol.GetNumAtoms() <= 3:
                has_3d = True
    results["conformer_3d"] = {
        "passed": has_3d,
        "message": "有效 3D 构象空间坐标" if has_3d else "仅包含 2D 平面坐标或缺失构象",
    }

    # 3. Internal clash check
    if mol.GetNumConformers() > 0 and has_3d:
        conf = mol.GetConformer()
        heavy = [i for i in range(mol.GetNumAtoms()) if mol.GetAtomWithIdx(i).GetSymbol() != "H"]
        bonded = {
            (min(b.GetBeginAtomIdx(), b.GetEndAtomIdx()), max(b.GetBeginAtomIdx(), b.GetEndAtomIdx()))
            for b in mol.GetBonds()
        }
        neighbors: dict[int, set[int]] = {}
        for b in mol.GetBonds():
            neighbors.setdefault(b.GetBeginAtomIdx(), set()).add(b.GetEndAtomIdx())
            neighbors.setdefault(b.GetEndAtomIdx(), set()).add(b.GetBeginAtomIdx())
        geminal = set()
        for nbs in neighbors.values():
            for a in nbs:
                for b in nbs:
                    if a < b:
                        geminal.add((a, b))

        pos = [conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]
        min_d = math.inf
        clash_pair = None
        for a_idx, i in enumerate(heavy):
            for j in heavy[a_idx + 1:]:
                if (i, j) in bonded:
                    continue
                thresh = 1.9 if (i, j) in geminal else 1.6
                d = math.dist((pos[i].x, pos[i].y, pos[i].z), (pos[j].x, pos[j].y, pos[j].z))
                if d < thresh and d < min_d:
                    min_d = d
                    clash_pair = (i, j)
        if clash_pair is not None:
            results["internal_clash"] = {
                "passed": False,
                "message": f"原子 {clash_pair[0]}-{clash_pair[1]} 碰撞 ({min_d:.2f} Å < 1.60 Å)",
            }
        else:
            results["internal_clash"] = {"passed": True, "message": "无配体内部原子严重空间位阻碰撞"}

    # 4. Protein clash check (when receptor is available)
    if rec_atoms and mol.GetNumConformers() > 0 and has_3d:
        conf = mol.GetConformer()
        heavy_pos = [
            conf.GetAtomPosition(i)
            for i in range(mol.GetNumAtoms())
            if mol.GetAtomWithIdx(i).GetSymbol() != "H"
        ]
        min_rec_d = math.inf
        for lp in heavy_pos:
            for ra in rec_atoms:
                d = math.dist((lp.x, lp.y, lp.z), ra["pos"])
                if d < min_rec_d:
                    min_rec_d = d
        if min_rec_d < 1.8:
            results["minimum_distance_to_protein"] = {
                "passed": False,
                "message": f"与受体原子碰撞 ({min_rec_d:.2f} Å < 1.80 Å)",
            }
        else:
            results["minimum_distance_to_protein"] = {
                "passed": True,
                "message": f"与受体间距合规 (最小距离 {min_rec_d:.2f} Å)",
            }

    # 5. Molecular formula check
    try:
        formula = rdMolDescriptors.CalcMolFormula(mol)
        results["molecular_formula"] = {"passed": True, "message": f"原子守恒 ({formula})"}
    except Exception:
        results["molecular_formula"] = {"passed": False, "message": "分子式计算失败"}

    return {"results": results}


def render_ascii_mode(
    target_path: Path,
    receptor_text: str,
    poses_list: list[str],
    options: dict[str, Any] | None = None,
) -> int:
    """Fallback terminal ASCII / interaction tree inspector."""
    console.print(make_header())

    first_pose_text = poses_list[0] if poses_list else ""
    affinity = _extract_affinity(target_path, first_pose_text)

    info_table = Table(box=None, show_header=False, expand=True)
    info_table.add_column("Key", style="bold #999999", width=14)
    info_table.add_column("Val", style="#FFFFFF")

    info_table.add_row("目标路径", str(target_path))
    info_table.add_row(
        "受体状态",
        f"已载入 ({len(receptor_text)} 字节)" if receptor_text else "[dim]未指定受体蛋白[/dim]",
    )
    info_table.add_row(
        "构象数量",
        f"{len(poses_list)} 个 Pose 构象" if poses_list else "[dim]无配体构象[/dim]",
    )
    if affinity is not None:
        info_table.add_row("结合亲和力", f"[bold #4EBA65]{affinity:.2f} kcal/mol[/bold #4EBA65]")
    else:
        info_table.add_row("结合亲和力", "[dim]未包含打分元数据 (Raw Structure)[/dim]")

    console.print(make_pane("终端构象概况 (ASCII Mode)", info_table))

    rec_atoms = _parse_receptor_atoms(receptor_text)
    lig_atoms, mol = _parse_ligand_atoms_and_mol(first_pose_text) if first_pose_text else ([], None)

    # Check for precalculated analysis in job directory
    interactions: list[dict[str, Any]] = []
    has_precalc = False
    if target_path.is_dir():
        p1 = target_path / "analysis" / "pose_1.json"
        if p1.is_file():
            try:
                data = json.loads(p1.read_text(encoding="utf-8"))
                interactions = data.get("interactions", [])
                has_precalc = True
            except Exception:
                pass

    if rec_atoms and lig_atoms:
        if not has_precalc:
            interactions = _calculate_dynamic_interactions(rec_atoms, lig_atoms)
        tree = render_interaction_tree(
            receptor_name="Receptor" if not target_path.is_dir() else target_path.name,
            ligand_name=target_path.stem if not target_path.is_dir() else "Pose_1",
            interactions=interactions,
            pose_index=1,
            affinity=affinity if affinity is not None else 0.0,
        )
        console.print(tree)
        console.print()
    elif mol is not None or lig_atoms:
        if mol is not None:
            console.print(_render_ligand_properties_table(mol, target_path.stem))
            console.print()
        else:
            console.print("[dim]已载入配体坐标，但未成功生成拓扑对象。[/dim]\n")
    elif rec_atoms:
        console.print(_render_receptor_summary_table(rec_atoms, target_path.stem))
        console.print("[#999999]当前仅载入受体蛋白结构，未提供配体构象，无相互作用力分析。[/]\n")
    else:
        console.print("[bold #FF6B80]未能识别或解析有效的受体或配体分子数据 (支持 PDB, PDBQT, SDF)[/bold #FF6B80]\n")

    qc_data = _run_rdkit_qc_checks(mol, rec_atoms if rec_atoms else None)
    console.print(render_posecheck_table(qc_data))
    return 0


def run_view(
    path: str | Path,
    port: int = 0,
    ascii_mode: bool = False,
    options: dict[str, Any] | None = None,
) -> int:
    """Runs the 3D viewer or terminal ASCII inspector.

    Returns:
        0: Normal exit on Enter or Ctrl+C, or successful ASCII rendering.
        1: Runtime or server error.
        2: Target structure or directory not found.
    """
    opts = options or {}
    target = Path(path).expanduser().resolve()
    settings = get_settings()

    receptor_text = ""
    poses_list: list[str] = []

    # 1. Resolve structures
    if target.is_dir():
        # Look for out.pdbqt, receptor.pdbqt, poses/*.sdf
        out_pdbqt = target / "out.pdbqt"
        rec_pdbqt = target / "receptor.pdbqt"
        rec_pdb = target / "receptor.pdb"
        if rec_pdbqt.is_file():
            receptor_text = rec_pdbqt.read_text(encoding="utf-8", errors="replace")
        elif rec_pdb.is_file():
            receptor_text = rec_pdb.read_text(encoding="utf-8", errors="replace")

        if out_pdbqt.is_file():
            poses_list = _split_models(out_pdbqt.read_text(encoding="utf-8", errors="replace"))
        else:
            sdf_files = sorted(
                target.glob("*.sdf"),
                key=lambda p: int(re.search(r"\d+", p.stem).group()) if re.search(r"\d+", p.stem) else 0,
            )
            if not sdf_files:
                sdf_files = sorted(
                    (target / "poses").glob("*.sdf"),
                    key=lambda p: int(re.search(r"\d+", p.stem).group()) if re.search(r"\d+", p.stem) else 0,
                )
            poses_list = [f.read_text(encoding="utf-8", errors="replace") for f in sdf_files]
    elif target.is_file():
        text = target.read_text(encoding="utf-8", errors="replace")
        if target.suffix == ".pdbqt":
            if "MODEL" in text:
                poses_list = _split_models(text)
            else:
                poses_list = [text]
            # Try to find sibling receptor
            sibling_rec = target.parent / "receptor.pdbqt"
            if sibling_rec.is_file():
                receptor_text = sibling_rec.read_text(encoding="utf-8", errors="replace")
        elif target.suffix == ".pdb":
            receptor_text = text
        elif target.suffix in (".sdf", ".mol", ".mol2"):
            poses_list = [text]
            # Try to find sibling receptor
            for cand in ["receptor.pdbqt", "receptor.pdb", "protein.pdb"]:
                sibling_rec = target.parent / cand
                if sibling_rec.is_file():
                    receptor_text = sibling_rec.read_text(encoding="utf-8", errors="replace")
                    break
    else:
        # Check ProjectStore for job ID
        store = ProjectStore(settings)
        job = store.find_job(str(path))
        if job is not None:
            job_dir = store.job_dir(job.project_id, job.id)
            return run_view(job_dir, port=port, ascii_mode=ascii_mode, options=options)

        console.print(f"\n[bold {COLOR_ERROR}]文件或路径不存在:[/bold {COLOR_ERROR}] '{path}'\n")
        return 2

    # 2. Terminal ASCII fallback
    if ascii_mode or opts.get("ascii"):
        return render_ascii_mode(target, receptor_text, poses_list, options=opts)

    # 3. Dynamic HTTP Server Setup
    actual_port = find_available_port(preferred_port=8080 if port == 0 else port)
    html_page = build_viewer_html(
        receptor_text=receptor_text,
        poses_list=poses_list,
        title=f"DockStudio - {target.name}",
        subtitle=target.name,
    )
    handler_class = create_viewer_handler(html_page)

    try:
        server = ReusableThreadingHTTPServer(("127.0.0.1", actual_port), handler_class)
    except OSError as exc:
        console.print(f"[bold {COLOR_ERROR}]无法绑定端口 {actual_port}:[/bold {COLOR_ERROR}] {exc}")
        return 1

    server_thread = threading.Thread(target=server.serve_forever, daemon=True)
    server_thread.start()

    # Determine SSH command
    try:
        user = getpass.getuser()
    except Exception:
        user = "username"
    host = socket.gethostname() or "remote-host"
    ssh_tunnel_cmd = f"ssh -N -L {actual_port}:localhost:{actual_port} {user}@{host}"

    # Render Announcement Pane
    announcement = Text()
    announcement.append("• 本地访问地址: ", style="#999999")
    announcement.append(f"http://127.0.0.1:{actual_port}\n\n", style="bold #4EBA65")

    announcement.append("• SSH 隧道转发命令 (在本地笔记本终端复制运行):\n", style="bold #93A5FF")
    announcement.append(f"  {ssh_tunnel_cmd}\n\n", style="bold #FFFFFF")

    announcement.append("• 转发建立后本地浏览器打开: ", style="#999999")
    announcement.append(f"http://localhost:{actual_port}\n\n", style="bold #48968C")

    announcement.append("• 运行状态: ", style="#999999")
    announcement.append(f"已就绪 ({len(poses_list)} 个构象) · 即用即焚微型服务\n", style="#CCCCCC")

    console.print(make_header())
    console.print(make_pane(
        "DockStudio 微型 3D 分子构象预览服务",
        announcement,
        subtitle="按 [Enter] 或 [Ctrl+C] 退出并释放端口",
        title_style="bold #D77757",
    ))

    # Wait for exit trigger
    try:
        input()
    except (KeyboardInterrupt, EOFError):
        pass
    finally:
        server.shutdown()
        server.server_close()
        server_thread.join(timeout=1.0)
        console.print(f"[dim]✓ 3D 预览服务已安全终止，端口 {actual_port} 已释放[/dim]\n")

    return 0
