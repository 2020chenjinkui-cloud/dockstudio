"""DockStudio backend package."""

# Linux manylinux wheels 的 ABI 冲突防护（实测：ubuntu CI / docker-slim 会段错误，
# Windows wheels 无此问题）：
# 1. openbabel 必须先于 rdkit 加载会段错误 → 包入口强制先加载 rdkit；
# 2. rdkit 部分子模块（MolStandardize/InChI/posebusters 链路）若在使用现场才
#    惰性导入，与已加载的 openbabel 库冲突同样段错误（exit 139）——因此把分析链路
#    会用到的 rdkit 子模块与 posebusters 也在此处（openbabel 之前）一并预热。
# 注意顺序：本文件的 imports 必须保持 rdkit 系在包内任何 openbabel 引入之前。
from rdkit import Chem as _Chem  # noqa: F401
from rdkit.Chem import AllChem as _AllChem  # noqa: F401


def _optional_preload() -> None:
    """分析链路的惰性导入预热（缺包无害，仅最佳努力）。"""
    import importlib

    for mod in (
        "rdkit.ForceField",
        "rdkit.ForceField.rdForceFieldHelpers",
        "rdkit.Chem.MolStandardize",
        "rdkit.Chem.SaltRemover",
        "rdkit.Chem.rdinchi",
        "posebusters",  # 其内部惰性 import MolStandardize，提前整包预热
    ):
        try:
            importlib.import_module(mod)
        except Exception:  # noqa: BLE001
            continue


_optional_preload()

__version__ = "2.0.3"
