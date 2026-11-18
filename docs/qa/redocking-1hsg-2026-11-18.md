# 1HSG 共晶配体自对接（redocking）RMSD 实测报告

| 项 | 值 |
|---|---|
| 任务 | dockstudio-qa-round1 / **t3** — 1HSG 真实 redocking RMSD 实测并对齐 D9 判据 ≤ 2.0 Å |
| 执行 | redock-verifier（attempt 1，2026-11-18） |
| 实测窗口 | 16:37–17:23 CST（5 批，共 14 次真 vina 对接） |
| 被测代码 | `main` @ `c0f33cc7dbea2928d9840a18f4ab11c66d8bc317`（DockStudio 2.0.5） |
| 测试工作树 | `/home/kui/wt/redock`（分支 `qa/redock-1hsg`；**未修改 `server/dockstudio/**`**） |
| 测试跑道 | `kui@100.98.223.80`（Tesla P4 8 GiB；Vina 为纯 CPU 引擎，GPU 未参与） |
| 产物 | `/home/kui/qa-runs/redock-1hsg/`（`results.json`＝10 个 run 全量、`api-check.json`、`ligand-prep-sweep.json`、各 run 的 `work/`） |

---

## 0. 结论（先看这里）

**D9 数值判定：PASS。** 用**共晶配体真实立体异构体**在 P4 上走平台代码路径做 1HSG 自对接：

| 口径（函数） | 数值 | 判据 |
|---|---|---|
| **top1 pose RMSD** — `analysis/crystal.py::heavy_atom_rmsd`（= `rdMolAlign.GetBestRMS`，对称校正 **+ 叠合**） | **0.380 Å** | ≤2.0 ✅ |
| **best-of-N RMSD（N=9）** — 同上口径取最小 | **0.380 Å** | ≤2.0 ✅ |
| **symmetry-corrected（`rdMolAlign.CalcRMS`，原位不叠合）** — D9「原位点」语义 | **0.454 Å** | ≤2.0 ✅ |

`exhaustiveness=32` 复测：0.553 / 0.428 / 0.598 Å，同样 PASS。**当前 2.0.5 + vina 1.2.7 在 1HSG 上没有 redocking 精度问题。**

**但 D9 的"证据有效性"不通过**，三处结构性缺陷会让不达标的对接显示为绿色：

1. **度量选错了函数。** `heavy_atom_rmsd` 底层 `GetBestRMS` 先做最优刚体叠合再算 RMSD。实测把正确姿态整体平移 +10 Å，平台值仍报 **0.380 Å**（`CalcRMS` 则 0.454 → **9.903 Å**）。本次 27 个姿态里 **6 个**在平台口径下 ≤2 Å（最低 0.428 Å）而原位偏差 >10 Å。而 D9 原文要求的是"**在原位点**重新对接的 RMSD"。
2. **用来 redock 的配体不是共晶配体。** 产品示例 `/api/demo/example`（也是仓库测试、也是两次历史 E2E）用的 `demo_data/indinavir_smiles.txt` 是**无立体** SMILES，实际落地的是 MK1 的**非对映体**（`CBVCZFGXHXORBI-RXFGFDIZSA-N` vs 晶体 `…-PXQQMZJSSA-N`，2 个手性中心相反），却照样报 1.075 Å"通过"。更严重的是**换一个 ETKDG 随机种子就会换一个分子**（见 §6.3：只测了 4 个种子，就给出 4 个**互不相同**的 InChIKey，其中只有 1 个恰好等于晶体 MK1）。
3. **证据不可复现。** `DockParams.seed` 默认 `None`（`models.py:31`）→ demo/E2E 的 Vina 搜索用随机种子；两次历史记录都没有记录 seed / vina 版本 / exhaustiveness / 度量函数。

**历史 4.226 Å（Windows 源码）/ 2.104 Å（macOS 安装包）的差异来源已定位到机制**（§6）：同机同版本下，改 **Vina 种子**只能把 RMSD 动到 1.08–1.22 Å（6 次运行 / 4 种种子设置），改 **配体准备环境** 则能把 RMSD 直接推到 **2.34 Å（叠合）/ 4.15 Å（原位）**——因为无立体 SMILES 的 3D 落地结果取决于 RDKit/ETKDG 的随机指派，不同平台的 RDKit 版本会给出不同的非对映体。即：**历史超标的主因是"配体身份不确定"，不是引擎或参数。**

---

## 1. 环境与版本（复现锚点）

| 项 | 值 |
|---|---|
| OS / 内核 | Linux 7.0.0-38-generic x86_64, glibc 2.43 |
| CPU | nproc = 4 |
| Python | 3.11.15（`/home/kui/miniconda3/envs/ai-gpu/bin/python`） |
| DockStudio | 2.0.5，实际导入 `/home/kui/wt/redock/server/dockstudio/__init__.py`（工作树，与 main 同一 commit） |
| RDKit | 2026.03.6 |
| meeko | 0.8.0 |
| OpenBabel（pose 导出用） | 3.2.1（python 绑定） |
| numpy | 2.4.4 |
| Vina | **1.2.7**，`/home/kui/dockstudio/resources/bin/linux-x86_64/vina`<br>sha256 `f31f774f723bba7bbe6e9d1c47577020eea9a8da16424284c043d22593570644` |

> 工作树 `/home/kui/wt/redock/resources/bin/` 不存在（`resources/bin/` 被 `.gitignore`），所以运行时显式指定 `DOCKSTUDIO_VINA=…`；该二进制即主仓库 2.0.5 自带的那一份（sha256 已记录，见附录 A）。

---

## 2. 平台代码路径：RMSD 到底怎么算的

### 2.1 被测函数

`server/dockstudio/analysis/crystal.py:34-47`

```python
def heavy_atom_rmsd(crystal_sdf: Path, pose_sdf: Path) -> float | None:
    crystal = _read_sdf(crystal_sdf)                  # MolFromMolBlock(sanitize=False, removeHs=False)
    pose    = _read_sdf(pose_sdf)
    ref   = _strip_bond_types(Chem.RemoveHs(crystal))  # 键全改 SINGLE、清芳香（crystal.py:22）
    probe = _strip_bond_types(Chem.RemoveHs(pose))
    rmsd  = rdMolAlign.GetBestRMS(probe, ref)          # ← 关键
    return round(float(rmsd), 3)
```

`GetBestRMS` 同时做两件事：
1. 枚举参考分子图的对称自同构并取最小 → **symmetry-corrected = 是**；
2. **先做最优刚体叠合**再算 RMSD → **原位（in-place）偏差 = 否**。

RDKit 自己把两种口径分得很清楚（本机 `CalcRMS.__doc__` 原文）：

> `CalcRMS`: "In contrast to `getBestRMS`, the RMS is computed **'in place'**, i.e. probe molecules are not aligned to the reference ahead of the RMS calculation. This is useful, for example, to compute the RMSD between **docking poses and the co-crystallized ligand**."

即：**RDKit 明确推荐 `CalcRMS` 做 redocking 判据，平台用的是 `GetBestRMS`。**

### 2.2 调用链（与产品行为一一对应）

| 环节 | 位置 | 行为 |
|---|---|---|
| 共晶 SDF 生成 | `prep/receptor.py:296-344` | PDB `HETATM` 分组 → `MolFromPDBBlock(sanitize=False)` → `MolToMolBlock`。**无键级**（PDB 不带键级），只有坐标 + 按距离感知的连接性 |
| 姿态导出 | `analysis/export.py:33-66` `export_poses` | `out.pdbqt` 按 `MODEL/ENDMDL` 切块，**OpenBabel(pybel)** 写 `pose_<n>.sdf`（顺序 = vina 输出顺序 = 亲和力排序） |
| 挂到 Job | `core/queue.py:254-258` | `pose.file = pose_<n>.sdf` |
| HTTP 端点 | `api/app.py:918-938` `GET /api/jobs/{jid}/crystal_rmsd?pose_index=N` | `pose_paths[N-1]` → `heavy_atom_rmsd(crystal_sdf, pose)` |
| 前端徽标 | `apps/web/src/viewerLogic.ts:26-29` `crystalBadge` | `rmsd < 2 → 'pass'`（绿色"通过"） |
| 另一处实现 | `cli/commands/verify.py:302-317` `calculate_best_rmsd` | 同样 `GetBestRMS`；图不匹配时退化为"按原子序号的直接坐标 RMSD" |

### 2.3 本报告同时给出的三种口径

| 字段名（`results.json`） | 函数 | 对称校正 | 叠合 |
|---|---|---|---|
| `platform_heavy_atom_rmsd` | `rdMolAlign.GetBestRMS`（经平台函数） | 是 | **是** |
| `CalcRMS_sym_inplace` | `rdMolAlign.CalcRMS` | 是 | **否（原位）** |
| `raw_firstmatch_inplace` | `prb.GetSubstructMatch(ref)` 首匹配后直接按坐标算 | 否 | 否 |

三者都在 `RemoveHs` + "所有键改 SINGLE/去芳香"之后计算，与平台口径的原子集一致（见 `scripts/qa_redock_1hsg.py::rmsd_suite`）。

---

## 3. 输入、盒子与配体身份校验

### 3.1 输入文件

| 输入 | 路径 | 校验 |
|---|---|---|
| 受体 | `server/dockstudio/demo_data/1HSG.pdb`（`/api/demo/example` 用的那份） | md5 `e0a74a5320205868bd0ffcad7fddb8ed`，与 `server/tests/fixtures/1HSG.pdb`、`fixtures/golden/1HSG.pdb` **完全一致** |
| 共晶配体参考 | `prepare_receptor` 自动抽取 `MK1:B` | 45 个重原子；SDF 3865 B（与 macOS E2E 记录 `crystal_sdf_bytes: 3865` 一致） |
| 配体（示例 / 仓库测试） | `demo_data/indinavir_smiles.txt`（= `fixtures/indinavir_smiles.txt`） | md5 `cc6a8ccbf3e13163b2751ecf8877341e` |
| 配体（正确立体） | `server/tests/golden/targets.json` → `pdb_id=1HSG, ccd=MK1` 的 `smiles` | 仓库自带金标准口条 |

### 3.2 盒子（对应排查项"盒子是否覆盖共晶配体位点/质心偏差"→ **排除**）

| 项 | 值 | 来源 |
|---|---|---|
| 中心 | (13.073, 22.467, 5.557) | `PreparedReceptor.cocrystal_ligands[0].center`（MK1:B 重原子质心，`prep/receptor.py:296-316`） |
| 与晶体配体质心偏差 | **0.001 Å** | 独立重算 |
| 边长 | 22.5 × 22.5 × 22.5 Å | demo 端点 / `test_audit.py` 常量 |
| 晶体配体原子是否全在盒内 | **45/45** | 最大轴偏差 x 8.87 / y 5.11 / z 4.62 Å < 半边长 11.25 Å |

盒子与 macOS E2E 记录的 `box.center=(13.073,22.467,5.557), size=(22.5,22.5,22.5)` 逐位相同。

### 3.3 配体身份校验（关键发现之一）

"从 PDB 提共晶配体 → 用金标准模板补键级 → `AssignStereochemistryFrom3D` 从晶体坐标反推立体化学"：

| 分子 | 来源 | InChIKey | 骨架块 | 立体块 |
|---|---|---|---|---|
| **晶体 MK1:B（基准）** | 1HSG HETATM + 模板补键级 + 3D 反推立体 | `CBVCZFGXHXORBI-PXQQMZJSSA-N` | CBVCZFGXHXORBI | `PXQQMZJSSA` |
| 仓库金标准 SMILES | `golden/targets.json` | `CBVCZFGXHXORBI-PXQQMZJSSA-N` | 同 | 同（**与晶体一致**） |
| **产品示例 / 仓库测试 SMILES** | `demo_data/indinavir_smiles.txt` | `CBVCZFGXHXORBI-UHFFFAOYSA-N` | 同 | `UHFFFAOYSA`（**未指定**） |
| **示例 SMILES 实跑出来的 3D 分子** | Run A 全部 9 姿态坐标反推 | `CBVCZFGXHXORBI-RXFGFDIZSA-N` | 同 | `RXFGFDIZ`（**另一非对映体**） |

两处手性中心相反（同骨架 canonical SMILES）：

```
晶体 MK1 : CC(C)(C)NC(=O)[C@@H]1CN(Cc2cccnc2)CCN1C[C@@H](O)C[C@@H](Cc1ccccc1)C(=O)N[C@H]1c2ccccc2C[C@H]1O
示例对接 : CC(C)(C)NC(=O)[C@H]1CN(Cc2cccnc2)CCN1C[C@@H](O)C[C@H](Cc1ccccc1)C(=O)N[C@H]1c2ccccc2C[C@H]1O
                                  ^^                                            ^^
```

组成式 C36H47N5O4、45 重原子、11 可旋转键**全都一致**，所以元素/尺寸/可旋转键这类常规 QC 发现不了。`/api/demo/example` 的注释写"配体为 indinavir（SMILES 来自 PDB 共晶配体）"，与事实不符（只是同一组成式）。

### 3.4 由共晶 SDF 反建配体会静默变成另一个分子（关键发现之二）

把平台自己抽出的 `cocrystal_MK1_B.sdf` 直接喂给 `prepare_ligand_from_file`（"用共晶配体做真正的 self-redock"最直觉的做法），**没有任何报错/警告**，得到的是：

| | 正确 MK1 | 由共晶 SDF 反建 |
|---|---|---|
| 组成式 | C36H47N5O4 | **C36H69N5O4**（多 22 个 H） |
| 可旋转键 | 11 | **13** |
| InChIKey | `CBVCZFGXHXORBI-PXQQMZJSSA-N` | `LCBGCXTWZMACCR-ZTFFRWSKSA-N` |
| 实际结构 | 含 3 个芳香环 + 1 个吡啶、1 个酰胺 C=O | **芳香环全部饱和成环己烷/哌啶，C=O 变成 C–OH** |

原因：共晶 SDF 无键级（§2.2 第一行），`prep/ligand.py` 读入后 `sanitize=False → SanitizeMol` 能通过（全单键价键合法），ETKDG 按饱和环嵌入。实测该分子的对接结果直接 FAIL（§5.4：top1 3.641 Å、best 2.037 Å、原位 9.660 Å）。产品当前没有地方这么用，但它是"从共晶配体建 ligand"这条最自然路径上的静默陷阱。

### 3.5 配体准备参数

`prep/ligand.py::prepare_ligand_from_smiles` → `_prepare_core_mol`（脱盐 + 互变异构/电荷标准化）→ `prep/tcm.py::generate_conformers_and_min_energy`：

| 项 | 值 |
|---|---|
| 3D 生成 | `AllChem.ETKDGv3()`，`useMacrocycleTorsions=True`、`useSmallRingTorsions=True`，`randomSeed=42`（**硬编码**），`numConfs=10` |
| 力场 | MMFF94s，500 步最小化，取**最低能单构象**进入对接 |
| 观测 `min_energy_kcal` | 示例 SMILES 175.84 / 金标准 SMILES 184.59（不同立体异构体，能量不可比） |
| 可旋转键 | 11 |
| PDBQT 写出 | meeko 0.8.0 `MoleculePreparation` + `PDBQTWriterLegacy`，87 行 |
| 受体准备 | **meeko 路线成功**，`warnings=[]`（未退化到 OpenBabel 回退），`receptor.pdbqt` 1842 原子 |

---

## 4. 参数表（逐项可复现）

| 参数 | Run A | Run B | Run C | 来源 / 默认 |
|---|---|---|---|---|
| 配体 | 示例 SMILES（无立体，实为非对映体） | 金标准立体 SMILES | 金标准立体 SMILES | §3.1 |
| `box.center` | (13.073, 22.467, 5.557) | 同 | 同 | MK1:B 质心 |
| `box.size` | (22.5, 22.5, 22.5) | 同 | 同 | demo / 测试常量 |
| `exhaustiveness` | 8 | 8 | **32** | `DockParams` 默认 8 |
| `num_modes` | 9 | 9 | 9 | `DockParams` 默认 9 |
| `energy_range` | 3.0 | 3.0 | 3.0 | `DockParams` 默认 3.0 |
| `cpu` | 4 | 4 | 4 | 本机 nproc=4 |
| **`seed`** | **42（显式）** | **42** | **42** | 显式固定；默认是 `None`（见 §6） |
| `scoring` | vina | vina | vina | `DockParams` 默认 |
| 引擎 | vina 1.2.7 | 同 | 同 | §1 |
| wall time | 91.2 s | 99.6 s | 286.2 s | 4 核独占 |

Vina 实际命令行（`engines/vina.py:244` `_build_cmd`；cwd=run/work，输入已复制为 `receptor.pdbqt`/`ligand.pdbqt`）：

```
vina --receptor receptor.pdbqt --ligand ligand.pdbqt \
  --center_x 13.073 --center_y 22.467 --center_z 5.557 \
  --size_x 22.5 --size_y 22.5 --size_z 22.5 \
  --exhaustiveness 8 --num_modes 9 --energy_range 3.0 \
  --cpu 4 --scoring vina --out out.pdbqt --seed 42
```

---

## 5. 结果

### 5.1 三个报数（任务要求）

| 数值 | **Run B（正确立体，D9 判定用）** | Run A（产品示例配体） | Run C（exh=32） |
|---|---|---|---|
| top1 pose RMSD（`heavy_atom_rmsd` / GetBestRMS） | **0.380 Å** | 1.075 Å | 0.553 Å |
| best-of-N RMSD（N=9，同上口径取最小） | **0.380 Å** | 1.048 Å | 0.428 Å |
| symmetry-corrected（`rdMolAlign.CalcRMS`，原位） | **0.454 Å** | 1.104 Å | 0.598 Å |
| 平台判定（<2 Å） | ✅ | ✅ | ✅ |
| 原位判定（≤2 Å） | ✅ | ✅ | ✅ |

### 5.2 逐姿态明细（A / B / C；`plat`=GetBestRMS 叠合，`inpl`=CalcRMS 原位，`raw`=无对称原位，`cd`=质心距晶体 Å）

**Run A — 示例 SMILES（产品 demo 路径），`exhaustiveness=8`，`seed=42`**

| pose | affinity | plat | inpl | raw | cd |
|---|---|---|---|---|---|
| 1 | −10.511 | **1.075** | 1.104 | 1.465 | 0.179 |
| 2 | −10.493 | 1.928 | **11.321** | 11.353 | 1.129 |
| 3 | −10.386 | 1.048 | **10.693** | 10.726 | 0.814 |
| 4 | −10.286 | 2.516 | 2.995 | 3.036 | 1.371 |
| 5 | −9.964 | 2.005 | 2.249 | 2.429 | 0.960 |
| 6 | −9.524 | 2.183 | 10.899 | 10.937 | 1.030 |
| 7 | −9.356 | 3.847 | 4.428 | 4.436 | 0.576 |
| 8 | −9.118 | 2.583 | 11.255 | 11.269 | 5.250 |
| 9 | −9.093 | 4.073 | 7.645 | 7.673 | 5.844 |

**Run B — 金标准立体 SMILES，`exhaustiveness=8`，`seed=42`**

| pose | affinity | plat | inpl | raw | cd |
|---|---|---|---|---|---|
| 1 | −11.265 | **0.380** | 0.454 | 0.846 | 0.150 |
| 2 | −10.557 | 0.428 | **10.921** | 10.960 | 0.745 |
| 3 | −10.151 | 3.888 | 4.171 | 4.211 | 0.573 |
| 4 | −10.053 | 3.977 | 10.810 | 10.843 | 1.870 |
| 5 | −9.996 | 2.427 | 11.389 | 11.420 | 2.022 |
| 6 | −9.955 | 2.299 | 4.223 | 4.259 | 3.331 |
| 7 | −9.778 | **1.811** | **10.986** | 11.025 | 0.730 |
| 8 | −9.710 | 3.881 | 10.820 | 10.849 | 1.992 |
| 9 | −9.628 | 2.247 | 10.836 | 10.877 | 1.390 |

**Run C — 金标准立体 SMILES，`exhaustiveness=32`，`seed=42`**

| pose | affinity | plat | inpl | raw | cd |
|---|---|---|---|---|---|
| 1 | −11.304 | **0.553** | 0.598 | 0.930 | 0.106 |
| 2 | −11.226 | 0.883 | 1.006 | 1.311 | 0.410 |
| 3 | −10.585 | 0.428 | **10.921** | 10.960 | 0.745 |
| 4 | −10.230 | 3.824 | 10.368 | 10.398 | 0.551 |
| 5 | −10.201 | 3.922 | 4.211 | 4.280 | 0.697 |
| 6 | −10.176 | 3.399 | 3.814 | 3.905 | 1.158 |
| 7 | −10.081 | 3.977 | 10.810 | 10.843 | 1.870 |
| 8 | −10.023 | 2.427 | 11.389 | 11.420 | 2.022 |
| 9 | −9.988 | 1.994 | 11.012 | 11.049 | 1.414 |

**统计（9 姿态 / run）**

| run | plat <2 Å | inpl <2 Å | 度量盲区（plat<2 且 inpl>5） |
|---|---|---|---|
| A | 3/9 | 1/9 | pose2 (1.928 vs 11.321)、pose3 (1.048 vs 10.693) |
| B | 3/9 | 1/9 | pose2 (0.428 vs 10.921)、pose7 (1.811 vs 10.986) |
| C | 4/9 | 2/9 | pose3 (0.428 vs 10.921)、pose9 (1.994 vs 11.012) |

**度量不变性决定性实验**（Run B 的 `pose_1.sdf` 整体沿 x 平移 +10 Å，内部几何完全不变）：

| 文件 | `heavy_atom_rmsd`（平台） | `CalcRMS`（原位） |
|---|---|---|
| `pose_1.sdf` 原样 | 0.380 Å | 0.454 Å |
| 同姿态 +10 Å | **0.380 Å（不变）** | **9.903 Å** |

→ 平台口径度量的是"构象对不对"，**不能证明"还在不在原位点"**；D9 原文（`交付物清单与验收标准_虚拟筛选项目.md` D9 行）要的是"在原位点重新对接"。

### 5.3 姿态选择与导出顺序（排查项"是否取了错误文件"→ **排除**）

- `out.pdbqt` 的 `REMARK VINA RESULT`：MODEL 1 恒为 `0.000 0.000`（相对本 run 最优模），MODEL 2+ 均 >0 → **MODEL 1 = 最佳亲和力模**。
- `export_poses` 按 `MODEL/ENDMDL` 顺序写 `pose_1..N.sdf`；`api` 取 `poses[N-1]`。
- 逐姿态核对 `pose_1.sdf` 亲和力/坐标 = MODEL 1（Run B：−11.265），**无跨姿态/跨文件错配**。
- 观察（不影响结论）：`export_poses` 走 OpenBabel 时打印 `Failed to kekulize aromatic bonds in OBMol::PerceiveBondOrders` —— 姿态 SDF 的**键级**可能不可靠；RMSD 路径先 `_strip_bond_types` 故不受影响，但下游读姿态 SDF 键级的模块（interactions / posecheck）需另测。

### 5.4 补充对照：由共晶 SDF 反建的配体（§3.4）

| pose | affinity | plat | inpl |
|---|---|---|---|
| 1 | −9.336 | **3.641** | **9.660** |
| 2 | −9.160 | 2.037 | 2.899 |
| 3 | −9.006 | 2.756 | 9.692 |
| best-of-9 | — | **2.037** | **2.899** |

即：喂错分子（饱和类似物）会让 D9 **FAIL**，但平台只会给出一个"3.641 Å"的数字，不会告诉你分子错了。

---

## 6. 与历史 4.226 Å / 2.104 Å 的差异来源（逐项证据）

### 6.1 历史记录的条件（本地只读仓库留档）

| 记录 | 位置 | 数值 | 同一步的其它观测 |
|---|---|---|---|
| Windows 源码模式 v2.0.4 | `DockStudio_代码审查与应用测试报告.md:58` | `crystal_rmsd` reference `MK1:B`，**4.226 Å** | `n_poses=5`、top **−9.855**、`/api/demo/example` `warnings=[]`（meeko 正常）、`avg_dock_seconds=156.43` |
| macOS 安装包 v2.0.4 | `e2e_macos_installed_app_20steps.json:145` | `{"reference":"MK1:B","rmsd":**2.104**,"crystal_sdf_bytes":3865}` | `n_poses=5`、top **−10.15**、demo box = (13.073, 22.467, 5.557) / 22.5（与本机逐位相同） |

两次记录**都没有**：seed、vina 版本、exhaustiveness/num_modes、度量函数、叠合/原位口径。

### 6.2 本机把"参数类"变量逐个隔离（`exhaustiveness` / `seed` / 复现性）

| run | 配体 | exh | **Vina seed** | n | top1 aff | plat top1 | plat best | inpl top1 |
|---|---|---|---|---|---|---|---|---|
| **A** | 示例 SMILES | 8 | 42 | 9 | −10.511 | 1.075 | 1.048 | 1.104 |
| **A2**（A 的重复运行） | 示例 SMILES | 8 | 42 | 9 | **−10.511** | **1.075** | **1.048** | **1.104** |
| **D** | 示例 SMILES | **2** | 42 | 3 | −10.713 | 1.185 | 1.185 | 1.284 |
| **E1** | 示例 SMILES | 8 | **None（默认）** | 9 | −10.796 | 1.209 | 1.079 | 1.299 |
| **E2** | 示例 SMILES | 8 | **None（默认）** | 9 | −10.557 | 1.077 | 1.077 | 1.106 |
| **F** | 示例 SMILES | 8 | **7** | 9 | −10.769 | 1.220 | 1.097 | 1.303 |
| **G** | 示例 SMILES | 8 | **2026** | 9 | −10.758 | 1.192 | 1.091 | 1.285 |
| **B** | 金标准立体 | 8 | 42 | 9 | −11.265 | **0.380** | 0.380 | **0.454** |
| **C** | 金标准立体 | 32 | 42 | 9 | −11.304 | 0.553 | 0.428 | 0.598 |
| **H** | 由共晶 SDF 反建（§3.4） | 8 | 42 | 9 | −9.336 | 3.641 | 2.037 | 9.660 |

三条结论（都是实测，不是推断）：

1. **固定 seed 后 vina 1.2.7 完全确定**：A 与 A2 的 top1 亲和力、top1/best 两个口径的 RMSD **逐位相同**（1.075 / 1.048 / 1.104，−10.511），wall 91.2 vs 90.0 s。
2. **Vina 种子不是历史超标的原因**：6 次运行（4 种种子设置：42 / 7 / 2026 / None）的 top1 全在 **1.077–1.220 Å**，原位 1.106–1.303 Å，无一接近 2 Å，更不用说 4.23 Å。但 **E1 ≠ E2** 说明默认 `seed=None` 时结果确实逐次不同（top1 亲和力差 0.24 kcal/mol），**不可复现**这一点成立。
3. **搜索深度不是原因**：`exhaustiveness` 2 / 8 / 32 的 top1 分别是 1.185 / 1.075 / 0.553 Å，全部 ≤2 Å（所以 §8 不做"靠加大 exhaustiveness 过关"的建议）。

### 6.3 差异的真正来源：配体身份 + 准备环境（补充实验，ETKDG 准备种子扫描）

对**同一个无立体示例 SMILES**，只改 ETKDG 的准备随机种子（产品里硬编码 42），Vina 侧 seed 固定 42、exh=8：

| ETKDG seed | 实际落地的 3D 分子 InChIKey | 与晶体一致？ | min E | top1 aff | **plat top1** | **inpl top1** |
|---|---|---|---|---|---|---|
| **42（产品默认）** | `…-RXFGFDIZSA-N` | ✗ 非对映体 | 175.84 | −10.511 | 1.075 | 1.104 |
| **1** | `…-QMFKYOCCSA-N` | ✗ 另一非对映体 | 175.89 | **−9.744** | **2.338** | **4.152** |
| **7** | `…-PXQQMZJSSA-N` | ✓ 就是晶体 MK1 | 177.57 | −11.529 | 0.899 | 0.940 |
| **2026** | `…-UILQWOSNSA-N` | ✗ 又一非对映体 | 179.32 | −11.154 | 0.898 | 0.988 |

**四个准备种子（42 / 1 / 7 / 2026）给出四个互不相同的立体异构体**，其中只有 seed=7 那个恰好等于晶体 MK1；对接结果随之在 0.899–2.338 Å（叠合）/ 0.940–4.152 Å（原位）之间摆动 —— 而 Vina 种子固定不变。

这条把历史数字钉死了：
- Windows 的 **4.226 Å** 与配体准备环境相关的量级（seed=1 一档：top1 2.34 Å 叠合 / 4.15 Å 原位，top1 亲和力 −9.744 vs 历史 −9.855）**同量级**；两次历史 run 的能量窗内姿态数只有 5（本机同参数 9）、top1 亲和力差 0.36–0.66 kcal/mol，也说明那两次的搜索/起始构象与本机不同。
- macOS 的 **2.104 Å** 落在同一区间的中段。
- 因此：**历史超标 = "无立体 SMILES + 不同平台的 RDKit/ETKDG 版本 → 落到不同的非对映体/起始构象"**，再加上未固定 Vina 种子与未记录口径。不是"引擎算错了"，也不是"1HSG 这个靶点做不好"。

> 说明（诚实边界）：产品把 ETKDG 种子硬编码为 42，所以"准备种子"不是用户可调参数；本实验用它来**模拟"不同 RDKit/ETKDG 构建环境"**。同一 RDKit 版本下结果稳定（A == A2），跨版本是否稳定**无法在本机验证**（本机只有一个 RDKit 2026.03.6）。但"无立体输入 → 立体指派由采样决定"这一因果已经由 4 个种子 → 4 个不同分子直接证明。

### 6.4 口径（叠合 vs 原位）的贡献

历史两次只报平台口径（GetBestRMS，叠合）。叠合只会让数值**偏小**，所以口径不能解释"历史偏大"；但它是"4.226 / 2.104 都算超标 vs 本机 1.075 算通过"这类**结论不可比**的根因之一：若改用原位口径，4.226 只会更大，而本机的 1.075 也会变成 1.104——量级不变，但语义从"构象像不像"变成"还在不在原位"。

### 6.5 平台口径 = 实际产品行为的交叉验证

用 `fastapi.testclient` 走真实 HTTP（与 `server/tests/test_audit.py::test_crystal_ligand_and_rmsd` 同一请求体：exh=2 / num_modes=3 / cpu=4 / seed=42 / 示例 SMILES / MK1 质心 22.5 Å 盒）：

| 来源 | pose1 | pose2 | pose3 | 越界 404 |
|---|---|---|---|---|
| `GET /api/jobs/{jid}/crystal_rmsd` | **1.185** | 2.554 | 4.030 | ✅ |
| 本脚本模块级同参数（Run D） | **1.185** | — | — | — |

→ 模块级测量脚本与产品端点**逐位一致**（1.185 Å），本报告主表的口径即产品口径。

---

## 7. 结论与 D9 适用性

### 7.1 数值结论

| 判据 | 结果 |
|---|---|
| D9：symmetry-corrected RMSD ≤ 2.0 Å（top1，**平台口径** GetBestRMS） | **PASS — 0.380 Å** |
| D9：symmetry-corrected RMSD ≤ 2.0 Å（top1，**原位口径** `CalcRMS`） | **PASS — 0.454 Å** |
| D9：best-of-N ≤ 2.0 Å（N=9） | **PASS — 0.380 / 0.454 Å** |
| 历史 4.226 Å / 2.104 Å 是否可复现 | **否**（机制已定位，见 §6.3） |

### 7.2 证据有效性结论

**(a) 度量函数选错。** D9 要"原位点"，平台用 `GetBestRMS`（叠合后算）。证据是实验而非推理：同姿态平移 10 Å，平台值 0.380 Å **不变**，原位值 9.903 Å；27 个姿态中 6 个"平台 ≤2 Å 而原位 >10 Å"。任何"位点漂移但构象正确"的失败样本都会显示绿色通过。

**(b) 用来 redock 的分子不是共晶配体。** 示例/测试/demo 库的 indinavir SMILES 无立体信息，落地 3D 是 MK1 的非对映体（§3.3），且换准备环境就会换分子（§6.3）；它能报 1.075 Å 是因为该非对映体在口袋里仍能摆出接近骨架的构象。**D9 会在"分子都不对"时判 PASS。** 同时，用共晶 SDF 反建配体会静默变成饱和类似物（§3.4），直接 FAIL 而平台不解释原因。

**(c) 证据不可复现。** `DockParams.seed` 默认 `None`（`models.py:31`）→ E1 ≠ E2；两次历史记录缺 seed / vina 版本 / exhaustiveness / 度量函数，导致同版本跨机器给出 4.226 / 2.104 / 1.075 三个数而无法裁决。

**综合判定：** D9 作为"数值门槛" **PASS**；作为"证据门槛" **FAIL**（口径 + 身份 + 可复现三处缺口）。补上 §8 的 R1–R3 之后，D9 才有资格作为"盒子与参数是否有效"的客观证据。

---

## 8. 最小修复建议（写给后续实现任务；本轮**未改**任何 `server/dockstudio/**`）

按"改动量 / 收益"排序。

| # | 问题 | 最小改动 | 落点 | 验收方式 |
|---|---|---|---|---|
| **R1** | 度量对刚体位移不敏感（§5.2、§7.2a） | `heavy_atom_rmsd` 同时返回原位口径：`rmsd_inplace = rdMolAlign.CalcRMS(probe, ref)`，保留现有 `GetBestRMS` 作为"构象匹配"字段；`CrystalRmsd` 增加 `rmsd_inplace`（向后兼容，不改现有字段语义）；前端徽标改用 `rmsd_inplace` | `server/dockstudio/analysis/crystal.py:34-47`、`server/dockstudio/models.py`（`CrystalRmsd`）、`api/app.py:918-938`、`apps/web/src/viewerLogic.ts:26-29` | 平移 10 Å 的姿态必须判 FAIL；`test_crystal_ligand_and_rmsd` 增断言 |
| **R2** | 示例/测试/demo 库配体是晶体配体的**非对映体**，且随环境漂移（§3.3、§6.3） | 把 `demo_data/indinavir_smiles.txt` 与 `demo_data/demo_library.smi` 的 indinavir 条目换成 `server/tests/golden/targets.json` 里已存的立体 SMILES（`…-PXQQMZJSSA-N`）；`fixtures/indinavir_smiles.txt` 同理 | `demo_data/indinavir_smiles.txt`、`demo_data/demo_library.smi`、`server/tests/fixtures/indinavir_smiles.txt` | 新增单测：示例配体 InChIKey == 晶体 MK1 InChIKey；改后 1HSG redock top1 应 ≈0.38 Å |
| **R3** | D9 链路没有"配体↔共晶配体"身份校验 | 在 `/api/demo/example` 与 `crystal_rmsd` 响应里加身份校验（骨架 InChIKey 块 + 立体块 + 手性中心数），不一致返回 `identity_match=false` + warning | `api/app.py:810-880`、`analysis/crystal.py` | 用无立体 SMILES 跑示例 → 必须出现"立体不一致"告警；用 §3.4 的共晶 SDF 反建 → 必须告警"键级/组成不一致" |
| **R4** | 默认随机种子（§6.2） | demo 请求体显式带 `seed=42`；或 `DockParams.seed` 默认改 42，另加显式 `seed=None` 表示"真随机" | `server/dockstudio/models.py:24-39`、`api/app.py` demo/E2E 请求体 | 连续 3 次同参 demo 的 top1 亲和力与 RMSD 完全一致 |
| **R5** | CI 门禁太松，超标仍绿 | `test_crystal_ligand_and_rmsd` 的 `0.0 < rmsd < 5.0` 收紧为 `rmsd_inplace <= 2.0`（并按 best-of-N 报数）；`tests/golden/test_golden.py` 目前**只断言亲和力区间、完全不测 RMSD**，建议加 per-target redock RMSD 断言（首版可用宽松阈值建基线） | `server/tests/test_audit.py:180-230`、`server/tests/golden/test_golden.py` | 用 §3.4 的饱和类似物配体（top1 3.641 Å）或把姿态平移 10 Å 时，CI 必须红 |
| **R6** | 报告/审计包不带环境指纹 | job result / audit 包记录 `engine_version`、vina `sha256`、`seed`、`exhaustiveness`、`box`、RMSD 口径 | `engines/vina.py::dock`、`core/queue.py`、`analysis/audit.py` | audit.zip 能读到以上 6 项 |
| **R7**（低优先） | 共晶 SDF 无键级 → `prepare_ligand_from_file` 静默产出饱和类似物（§3.4） | 该路径对"无键级 SDF"给出显式拒绝或强制要求模板 `AssignBondOrdersFromTemplate`；或在 `cocrystal_*.sdf` 生成时就写正确键级（从 CCD/模板） | `prep/receptor.py:318-344`、`prep/ligand.py::prepare_ligand_from_file` | 由共晶 SDF 反建的配体必须与晶体配体 InChIKey 一致或报错 |

**明确不建议做的事**：不要靠调大 `exhaustiveness`（或缩盒子）去让数字 ≤2 Å —— 本次已证明 `exhaustiveness=8` 就足够（0.380 Å），问题不在搜索预算。

---

## 9. 复现步骤（完整命令）

```bash
# 0) 工作树（不动主工作树）
cd /home/kui/dockstudio
git worktree add /home/kui/wt/redock -b qa/redock-1hsg main
cd /home/kui/wt/redock && git rev-parse HEAD     # c0f33cc7dbea2928d9840a18f4ab11c66d8bc317

# 1) 本报告新增的脚本（只在 scripts/ 与 docs/，不碰 server/dockstudio/**）
#    scripts/qa_redock_1hsg.py          主测量：prep → dock → export → 三口径 RMSD + 身份校验 + 盒子校验
#    scripts/qa_redock_api_check.py     HTTP 端点交叉核对（TestClient，逐 pose 调 crystal_rmsd）
#    scripts/qa_redock_ligand_sweep.py  补充实验：仅改 ETKDG 准备种子，隔离"配体身份 vs Vina 种子"
cd /home/kui/wt/redock
VINA=/home/kui/dockstudio/resources/bin/linux-x86_64/vina   # 工作树无 resources/bin，必须显式给
PY=/home/kui/miniconda3/envs/ai-gpu/bin/python
OUT=/home/kui/qa-runs/redock-1hsg

# 2) 主测量（全程持重活锁；含 A/A2/B/C/D/E1/E2/F/G/H 共 10 个 run）
DOCKSTUDIO_VINA=$VINA flock /home/kui/.dockstudio-heavy.lock -c \
  "$PY scripts/qa_redock_1hsg.py --repo /home/kui/wt/redock --out $OUT --cpu 4"

# 3) HTTP 端点交叉核对（与仓库测试同一条请求体）
DOCKSTUDIO_VINA=$VINA flock /home/kui/.dockstudio-heavy.lock -c \
  "$PY scripts/qa_redock_api_check.py /home/kui/wt/redock $OUT"

# 4) 补充实验（配体准备种子扫描）
DOCKSTUDIO_VINA=$VINA flock /home/kui/.dockstudio-heavy.lock -c \
  "$PY scripts/qa_redock_ligand_sweep.py /home/kui/wt/redock $OUT"

# 5) 结果
cat $OUT/results.json            # 环境指纹 + 10 个 run 的逐姿态三口径 RMSD + 姿态 3D 立体反推
cat $OUT/results-batch{1,2,3}.json
cat $OUT/api-check.json          # GET /api/jobs/{jid}/crystal_rmsd 的逐 pose 返回
cat $OUT/ligand-prep-sweep.json  # ETKDG 准备种子扫描
```

**单次运行的最小命令**（只复现 D9 那三个数，约 2 分钟）：

```bash
cd /home/kui/wt/redock && DOCKSTUDIO_VINA=/home/kui/dockstudio/resources/bin/linux-x86_64/vina \
  flock /home/kui/.dockstudio-heavy.lock -c \
  "/home/kui/miniconda3/envs/ai-gpu/bin/python scripts/qa_redock_1hsg.py \
     --repo /home/kui/wt/redock --out /home/kui/qa-runs/redock-1hsg --cpu 4 --only B_golden_exh8"
```

产物目录 `/home/kui/qa-runs/redock-1hsg/`：

| 路径 | 内容 |
|---|---|
| `data/prepared/receptors/<id>/` | `receptor.pdbqt`(1842 原子)、`cocrystal_MK1_B.sdf`(3865 B)、`cocrystals.json` |
| `data/prepared/ligands/<id>/ligand.pdbqt` | 各配体变体（fixture / golden / 由共晶 SDF 反建 / prep-seed 扫描） |
| `<run_id>/work/` | `out.pdbqt`、`pose_1..9.sdf`、`receptor.pdbqt`、`ligand.pdbqt` |
| `results*.json` | 全部实测数据（vina sha256、逐姿态三口径、姿态 3D 立体反推、盒子校验） |
| `api-check.json` / `ligand-prep-sweep.json` | 端点交叉核对 / 准备种子扫描 |

---

## 10. 覆盖范围与残余不确定性

**已覆盖**

- 平台真实链路：`prepare_receptor` → `prepare_ligand_from_smiles` → `VinaEngine.dock` → `export_poses` → `analysis.crystal.heavy_atom_rmsd`（与 `core/queue.py` 的 job 流程一一对应），并用 `TestClient` 走真实 HTTP 端点交叉核对（§6.5，1.185 == 1.185）。
- 真 vina 1.2.7（非 mock）、真 meeko 0.8.0 受体/配体准备（**非** OpenBabel 回退），共 14 次真实对接。
- 排查清单四项全部有证据结论：盒子覆盖（排除，0.001 Å / 45-45）、配体身份（**命中**，非对映体）、局部弛豫（排除，全局搜索 + exh 扫描）、姿态文件选择（排除，MODEL 顺序与导出逐位一致）。

**未覆盖 / 残余不确定性**

1. **Windows/macOS 二进制与 RDKit 版本差异无法在本机钉死**：两次历史记录未写 vina 版本与 RDKit 版本，本机只有 Linux vina 1.2.7 + RDKit 2026.03.6。§6.3 已证明机制（准备环境决定是哪个分子），但"历史那次具体是哪个非对映体"没有留档可查。
2. **准备种子扫描只做了 4 个种子**，足以证明"身份不确定"，但不构成分布的统计描述。
3. **未做多引擎对照**（gnina / Uni-Dock 属另一条线）。
4. **GPU 未参与**：Vina 是纯 CPU 引擎；Tesla P4 对 vina 无意义（对 gnina 才有意义）。
5. **`export_poses` 的 OpenBabel 键级**未单独验证（本报告结论不依赖它，但 §5.3 已记录告警）。

---

## 附录 A：关键数字一览

| 指标 | 值 |
|---|---|
| 晶体配体 MK1:B | 45 重原子 / `CBVCZFGXHXORBI-PXQQMZJSSA-N` |
| 示例 SMILES | `CBVCZFGXHXORBI-UHFFFAOYSA-N`（立体未定） |
| 示例 SMILES 实际对接分子 | `CBVCZFGXHXORBI-RXFGFDIZSA-N`（非对映体） |
| 盒子中心偏差 / 盒内原子 | 0.001 Å / 45-45 |
| vina 版本 / sha256 | 1.2.7 / `f31f774f723bba7bbe6e9d1c47577020eea9a8da16424284c043d22593570644` |
| **Run B（金标准配体, exh 8, seed 42）top1 / best / 原位** | **0.380 / 0.380 / 0.454 Å** |
| Run A（示例配体, exh 8, seed 42）top1 / best / 原位 | 1.075 / 1.048 / 1.104 Å |
| Run C（金标准配体, exh 32, seed 42）top1 / best / 原位 | 0.553 / 0.428 / 0.598 Å |
| Run H（共晶 SDF 反建的饱和类似物）top1 / best / 原位 | 3.641 / 2.037 / 9.660 Å |
| A vs A2（同 seed 重复） | 逐位相同 → 固定 seed 后确定 |
| 度量不变性（+10 Å 平移） | 平台 0.380 Å 不变；原位 0.454 → 9.903 Å |
| "平台 ≤2 Å 但原位 >5 Å" 的姿态 | 27 个中 6 个 |
| ETKDG 准备种子 42/1/7/2026 → 3D 分子 | 4 个种子给出 **4 个互不相同的 InChIKey**（仅 seed=7 恰为晶体 MK1）；plat top1 1.075 / 2.338 / 0.899 / 0.898 Å |
