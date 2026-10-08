"""DockStudio backend package."""

# Linux manylinux wheels 的 ABI 冲突：openbabel 先于 rdkit 导入时，随后再导入
# rdkit.Chem 会段错误（SIGSEGV，exit 139）。在包入口强制先加载 rdkit，
# 使任何 dockstudio.* 导入路径都保持安全顺序。Windows wheels 无此问题。
from rdkit import Chem as _Chem  # noqa: F401

__version__ = "0.1.0"
