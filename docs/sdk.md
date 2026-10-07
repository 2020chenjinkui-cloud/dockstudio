# DockStudio SDK（直接用 import 的调用方式）

SDK 就是 `dockstudio` Python 包本身 —— 各核心模块可以直接 import 使用，无需经过 HTTP API。

```bash
cd server && py -3.12 -m pip install -e .
```

## 最小示例：准备 + 对接 + 导出

```python
from pathlib import Path
from dockstudio.config import Settings
from dockstudio.prep.receptor import prepare_receptor
from dockstudio.prep.ligand import prepare_ligand_from_smiles
from dockstudio.engines.vina import VinaEngine
from dockstudio.analysis.export import export_poses
from dockstudio.models import Box, DockParams

settings = Settings(data_dir=Path("./my_data"))
receptor = prepare_receptor(Path("1HSG.pdb").read_bytes(), "1HSG.pdb", settings.prepared_dir / "receptors")
ligand = prepare_ligand_from_smiles("CCO", "ethanol", settings.prepared_dir / "ligands")

engine = VinaEngine()
box = Box(center=receptor.cocrystal_ligands[0].center, size=(22.5, 22.5, 22.5))
result = engine.dock(Path(receptor.pdbqt_path), Path(ligand.pdbqt_path), box, DockParams(), Path("./run"))
poses = export_poses(Path("./run/out.pdbqt"), Path("./run"))
print(result.poses[0].affinity)
```

## 相互作用分析与口袋预测

```python
from dockstudio.analysis.interactions import analyze_pose
from dockstudio.analysis.pockets import detect_pockets

analysis = analyze_pose(Path("receptor_protein.pdb"), Path("run/pose_1.sdf"), Path("run"))
pockets = detect_pockets(Path("receptor_protein.pdb"), top_n=3)
```

## 命令行入口（CLI）

```bash
py -3.12 -m dockstudio.cli serve      # 等价于 dockstudio:serve 入口，启动 uvicorn
py -3.12 -m dockstudio.cli engines    # 列出 engine manifest 与可用性
py -3.12 -m dockstudio.cli demo       # 用内置 1HSG fixture 跑完整流程并打印 top affinity
```

## 错误类型

所有模块抛出 `dockstudio.models` 里的异常：`PreparationError`、`DockingError`、`EngineUnavailableError`（均为 `DockStudioError` 子类，消息为用户可读文本）。
