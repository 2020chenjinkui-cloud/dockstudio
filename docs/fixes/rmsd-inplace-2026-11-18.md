# D9 度量口径修正报告：原位 RMSD（CalcRMS）成为判据，叠合 RMSD 降为辅助

- **日期**：2026-11-19
- **任务**：t19 —— 修正 D9 度量口径：新增原位 RMSD（CalcRMS）并区分叠合口径
- **分支**：`fix/rmsd-inplace` @ `/home/kui/wt/rmsd`（基线 `main` = `5f43063`）
- **未合 main、未 push**（按约定由 integrator 串行处理）
- **前置证据**：`docs/qa/redocking-1hsg-2026-11-18.md`（t3 的 1HSG 真机实测）

---

## 0. 摘要

| 项 | 修前 | 修后 |
|---|---|---|
| `heavy_atom_rmsd` 的度量 | `rdMolAlign.GetBestRMS`（**先叠合**再算距离） | 保留，但明确标注为**辅助口径**；不再充当 D9 结论 |
| D9 判据函数 | 无（用叠合口径冒充） | `rdMolAlign.CalcRMS`（**原位**，不叠合），导出为 `inplace_heavy_atom_rmsd` / `CrystalPoseRmsd.inplace_rmsd` |
| **正确姿态平移 +10 Å 后** | **0.0 Å（误判 PASS）** | **原位 9.903 Å → FAIL**；叠合仍 0.380 Å（辅助，不作判据） |
| `GET /api/jobs/{jid}/crystal_rmsd` | 只返回混淆的 `rmsd` | `rmsd`（原位，D9）+ `metric` + `rmsd_superposed` + `metric_superposed` + `threshold_a` + `d9_status` |
| HTML 报告 | 无共晶对照 | 新增 `redocking` 段：逐姿态双口径 + D9 判定 + best-of-N（按原位口径取最优） |
| audit.zip | 无共晶对照 | `report.md` 新增 D9 段；`snapshot.json` 新增机器可读 `cocrystal_rmsd` 块 |
| CLI `dockstudio audit` | 无共晶对照 | 同样打印 D9 段 |

**1HSG 真机复测（真 vina 1.2.7，P4）：原位 0.454 Å（PASS）/ 叠合 0.380 Å**，与 t3 预期一致；9 个姿态里只有 top1 通过 D9，其余 8 个原位 4.2–11.4 Å 全部 FAIL，而**其中 4 个的叠合值仍 < 4 Å**——这正是修前会被误判为"接近通过"的姿态。

---

## 1. 为什么要改：两个口径回答的是两个不同问题

D9 判据原文（`交付物清单与验收标准_虚拟筛选项目.md`）：

> **D9** | Redocking 对照 | … 共晶配体在**原位点**重新对接的 RMSD（目标 ≤ 2.0 Å）；**这是判断"盒子与参数是否有效"的唯一客观证据**

而修前 `analysis/crystal.py:34` 的 `heavy_atom_rmsd` 是：

```python
rmsd = rdMolAlign.GetBestRMS(probe, ref)     # 先最优刚体叠合，再算距离
```

`GetBestRMS` 在计算前会把 probe 叠合到 ref 上，于是**平移不变**：姿态漂出结合口袋，数字仍然是 0。用这样一个对刚体位移完全不敏感的度量去回答"盒子与参数是否有效"，逻辑上不成立。

RDKit 自己把两种口径分得很清楚（本机 `CalcRMS.__doc__` 原文）：

> `CalcRMS`: "In contrast to `getBestRMS`, the RMS is computed **'in place'**, i.e. probe molecules are not aligned to the reference ahead of the RMS calculation. This is useful, for example, to compute the RMSD between **docking poses and the co-crystallized ligand**."

### 1.1 两个口径的定义与公式

设参考（共晶配体）重原子坐标 $R = \{r_i\}$，姿态重原子坐标 $P = \{p_i\}$（均先 `RemoveHs` 并把所有键改 SINGLE，确保拓扑匹配不受键级感知差异影响）。令 $A$ 为参考分子图的自同构群（RDKit 内部枚举，含 `symmetrizeConjugatedTerminalGroups`）。

**原位 RMSD（D9 判据，`CalcRMS`）**

$$
\mathrm{RMSD_{inplace}} = \min_{a \in A} \sqrt{\frac{1}{N}\sum_{i=1}^{N}\lVert p_{a(i)} - r_i \rVert^2}
$$

$P$ 的坐标**一动不动**，只允许原子重排（对称校正）。这就是"还在不在原位点"。

**叠合 RMSD（辅助，`GetBestRMS`）**

$$
\mathrm{RMSD_{superposed}} = \min_{a \in A}\ \min_{T \in SE(3)} \sqrt{\frac{1}{N}\sum_{i=1}^{N}\lVert T\,p_{a(i)} - r_i \rVert^2}
$$

多了一个刚体自由度 $T$，于是对平移/旋转完全免疫——衡量的是"构象像不像"，不是"放对没放对"。

两者都做对称校正，所以**差值恰好等于"叠合"这一步带来的信息**：$\mathrm{RMSD_{superposed}} \ll \mathrm{RMSD_{inplace}}$ 就意味着"形状对、位置错"。

### 1.2 修前的实测后果（t3 数据，27 个姿态）

| 现象 | 数量 |
|---|---|
| 叠合口径 ≤ 2 Å（平台会判"通过"）而原位 > 10 Å | **6 / 27** |
| 最低的误判值 | **0.428 Å**（对应原位 10.921 Å） |

---

## 2. 实现

### 2.1 `analysis/crystal.py`：单一入口 + 明确命名

```python
INPLACE_METRIC   = "inplace_calc_rms"      # D9 判据
SUPERPOSED_METRIC = "superposed_best_rms"  # 辅助
CRITERION_METRIC = INPLACE_METRIC          # 明确"谁是判据"
AUXILIARY_METRIC = SUPERPOSED_METRIC
D9_THRESHOLD_A   = 2.0
```

| 导出物 | 作用 |
|---|---|
| `crystal_pose_rmsd(crystal_sdf, pose_sdf) -> CrystalPoseRmsd` | **推荐入口**。一次算出两个口径，各自独立捕获异常 |
| `CrystalPoseRmsd`（frozen dataclass） | `.inplace_rmsd` / `.superposed_rmsd` / `.metric` / `.auxiliary_metric` / `.rmsd`（=原位）/ `.d9_pass` / `.d9_status()` / `.as_dict()` |
| `inplace_heavy_atom_rmsd(...)` | 只要 D9 那一个数 |
| `superposed_heavy_atom_rmsd(...)` | 只要辅助那一个数 |
| `heavy_atom_rmsd(...)` | **保留兼容**，但 docstring 明确写"after best-fit superimposition … Deprecated as a D9 verdict"，且现在返回叠合值并提示新代码改用 `crystal_pose_rmsd` |
| `cocrystal_reference(receptor_dir)` | API / 报告 / 审计**共用**的共晶参考定位器（读 `cocrystals.json`） |
| `pose_rmsd_series(...)` / `best_pose_rmsd(...)` | 批量；`best_pose_rmsd` **按原位口径**取最优（按叠合取最优会挑出"形状好但跑出盒子"的姿态） |

两个口径各自 `try/except`，一方失败不影响另一方，且失败原因 `log.warning` 留痕（不静默降级）。

### 2.2 消费点：谁展示、谁进审计

| 消费点 | 拿到什么 | 位置 |
|---|---|---|
| `GET /api/jobs/{jid}/crystal_rmsd` | `rmsd`（原位）+ `metric` + `rmsd_superposed` + `metric_superposed` + `threshold_a` + `d9_status` + `reference` | `api/app.py` 端点 + 同文件内的 `CrystalRmsdResponse` |
| HTML 报告（`redocking` 段，默认开启） | 逐姿态：原位值 / D9 判定 / 叠合值（灰字标注"辅助"）+ best-of-N + 两个口径的全名与函数名 | `analysis/report.py:_redocking_html` |
| audit.zip `report.md` | 同上（Markdown 表）+ 明确写"叠合对刚体位移不敏感，不可作为 D9 结论" | `analysis/audit.py:build_audit_markdown` |
| audit.zip `snapshot.json` | 机器可读 `cocrystal_rmsd` 块：`criterion_metric` / `auxiliary_metric` / `threshold_a` / 每姿态 `as_dict()` | `analysis/audit.py:build_audit_zip` |
| CLI `dockstudio audit <job>` | 与 report.md 同一段（`build_audit_markdown` 是共享入口，CLI 调用它即自动获得 D9 段，无需改 CLI 文件） | `analysis/audit.py:build_audit_markdown` |

**字段命名规则（防混淆）**：凡是判据一律叫 `rmsd` / `inplace_rmsd`，辅助一律带 `superposed` 字样；每个数值旁边都能找到产生它的函数标识符，没有任何地方让两个口径共用同一个字段名。

> 兼容性：`CrystalRmsd.rmsd` 的字段名不变，但**含义从"叠合"变成"原位"**——这正是本次修复的目的（旧值是错的判据）。新增字段全部有默认值，前端现有的 `crystal?.rmsd` 读取路径不需要改就能拿到正确的 D9 数字（前端不在本任务 inScope）。

### 2.3 报告/审计性能

两个口径合计 **0.001 s/姿态**（9 姿态 0.01 s），报告与审计都可忽略不计。

---

## 3. 回归用例：先红后绿

### 3.1 红色证据（改实现之前）

用同一分子平移 +10 Å 的 SDF 对（苯乙酸 C8H8O2，10 重原子，ETKDGv3/MMFF94s 生成后整体 +10 Å）：

```
BEFORE: current heavy_atom_rmsd(crystal, pose+10A) = 0.0
BEFORE: same pose, untranslated                   = 0.0
```

即**平移 10 Å 完全不可见**。测试 `test_translated_pose_inplace_is_ten_angstrom` 当时报 `ImportError: crystal_pose_rmsd`，其余 5 条同因失败/报错。

### 3.2 绿色证据（改实现之后）

```
translated +10A -> {'rmsd': 10.0, 'metric': 'inplace_calc_rms',
                    'auxiliary_rmsd': 0.0, 'auxiliary_metric': 'superposed_best_rms',
                    'threshold_a': 2.0, 'd9_status': 'FAIL'}
identical       -> {..., 'rmsd': 0.0, ..., 'd9_status': 'PASS'}
```

`server/tests/test_crystal_rmsd_metric.py`（6 条，全绿）：

| 用例 | 断言 |
|---|---|
| `test_translated_pose_inplace_is_ten_angstrom` | 原位 9.5–10.5 Å |
| `test_translated_pose_superposed_is_near_zero` | 叠合 < 0.5 Å |
| `test_identical_pose_scores_zero_on_both` | 两者 < 0.05 Å |
| `test_metric_names_are_distinguishable` | `CRITERION_METRIC == INPLACE_METRIC != SUPERPOSED_METRIC`；`.rmsd is .inplace_rmsd` |
| `test_d9_pass_uses_inplace_metric` | 平移姿态 `d9_pass is False` |
| `test_heavy_atom_rmsd_docstring_states_its_metric` | 旧入口 docstring 声明自己是叠合口径 |

`server/tests/test_crystal_d9_wiring.py`（10 条，全绿）覆盖消费点：报告双口径+名称、best-of-N 按原位取、无共晶降级、API 双字段（含 null 场景仍带口径名）、audit.zip 的 report.md 与 snapshot.json、`pose_rmsd_series`/`best_pose_rmsd` 的"好形状但跑偏"姿态不得夺冠、`cocrystal_reference` 的降级路径。

---

## 4. 真机 1HSG 实测（P4，真 vina 1.2.7）

### 4.1 环境与参数（全部可复现）

| 项 | 值 |
|---|---|
| 机器 | `kui@100.98.223.80`，Tesla P4 8 GiB（vina 为纯 CPU 引擎，GPU 未参与），nproc=4 |
| OS / Python | Linux 7.0.0-38-generic x86_64 / CPython 3.11.15 |
| DockStudio | 2.0.5，`fix/rmsd-inplace` @ `5f43063f75e84da91fb10a18604059040eca8aa8`，实际导入 `/home/kui/wt/rmsd/server/dockstudio/__init__.py` |
| RDKit / meeko | 2026.03.6 / 0.8.0 |
| vina | **1.2.7**，`/home/kui/dockstudio/resources/bin/linux-x86_64/vina`<br>sha256 `f31f774f723bba7bbe6e9d1c47577020eea9a8da16424284c043d22593570644` |
| 受体 | `server/tests/fixtures/1HSG.pdb`（与 `demo_data/1HSG.pdb` md5 相同），`receptor.pdbqt` 1842 原子，`warnings=[]`（meeko 路线，无 OpenBabel 回退） |
| **盒子来源** | `PreparedReceptor.cocrystal_ligands[0].center` — MK1:B 45 个重原子的质心（`prep/receptor.py:296-316`） |
| **盒子中心** | **(13.073, 22.467, 5.557)** Å（与晶体配体独立重算质心的偏差 0.001 Å） |
| **盒子边长** | **(22.5, 22.5, 22.5)** Å（MK1 全部 45 个重原子均在盒内：max 轴向偏差 x 8.87 / y 5.11 / z 4.62 Å < 11.25 Å 半边长） |
| **exhaustiveness** | **8**（`DockParams` 默认） |
| **num_modes** | **9** |
| **seed** | **42**（显式固定；`DockParams.seed` 默认仍是 `None`，属另一任务） |
| energy_range / cpu / scoring | 3.0 / 4 / `vina` |
| 配体 SMILES | `server/tests/golden/targets.json` 的 1HSG/MK1 立体 SMILES（**不是** demo 的无立体 SMILES，见 §5） |
| 配体身份校验 | 准备后 InChIKey `CBVCZFGXHXORBI-PXQQMZJSSA-N` == 由晶体坐标（模板补键级 + `AssignStereochemistryFrom3D`）反推的 `CBVCZFGXHXORBI-PXQQMZJSSA-N` **一致** |
| 构象生成 | ETKDGv3 `numConfs=10` `randomSeed=42` + MMFF94s 500 步，取最低能单构象（min E = 184.593 kcal/mol）；11 可旋转键 |
| 质子化 | meeko 0.8.0 `PDBQTWriterLegacy`（氢来自 ETKDG/MMFF） |
| wall time | 123.8 s |

实际命令行（`engines/vina.py::_build_cmd`，cwd = run/work，输入已复制为裸文件名）：

```
vina --receptor receptor.pdbqt --ligand ligand.pdbqt \
  --center_x 13.073 --center_y 22.467 --center_z 5.557 \
  --size_x 22.5 --size_y 22.5 --size_z 22.5 \
  --exhaustiveness 8 --num_modes 9 --energy_range 3.0 \
  --cpu 4 --scoring vina --out out.pdbqt --seed 42
```

复现命令：

```bash
# 测量（持重活锁；产物 /home/kui/qa-runs/d9-inplace/measure.json）
cd /home/kui/wt/rmsd
DOCKSTUDIO_VINA=/home/kui/dockstudio/resources/bin/linux-x86_64/vina \
  flock /home/kui/.dockstudio-heavy.lock -c \
  "/home/kui/miniconda3/envs/ai-gpu/bin/python docs/fixes/rmsd-inplace-measure_1hsg_d9.py /home/kui/wt/rmsd /home/kui/qa-runs/d9-inplace"
```

### 4.2 结果：逐姿态双口径

| pose | affinity (kcal/mol) | **原位 RMSD (Å)** | **D9** | 叠合 RMSD (Å，辅助) |
|---|---|---|---|---|
| **1** | −11.265 | **0.454** | ✅ **PASS** | 0.380 |
| 2 | −10.557 | 10.921 | ❌ FAIL | 0.428 |
| 3 | −10.151 | 4.171 | ❌ FAIL | 3.888 |
| 4 | −10.053 | 10.810 | ❌ FAIL | 3.977 |
| 5 | −9.996 | 11.389 | ❌ FAIL | 2.427 |
| 6 | −9.955 | 4.223 | ❌ FAIL | 2.299 |
| 7 | −9.778 | 10.986 | ❌ FAIL | 1.811 |
| 8 | −9.710 | 10.820 | ❌ FAIL | 3.881 |
| 9 | −9.628 | 10.836 | ❌ FAIL | 2.247 |

**汇总**

| 口径 | top1 | best-of-9（按该口径取最优） |
|---|---|---|
| **原位（D9 判据）** | **0.454 Å → PASS** | **0.454 Å（pose 1）→ PASS** |
| 叠合（辅助） | 0.380 Å | 0.380 Å（pose 1） |

与 t3 的预期逐位吻合（原位 0.454 / 叠合 0.380）。**注意 pose 2 与 pose 7**：叠合只有 0.428 / 1.811 Å（若按修前口径会被认为"接近通过/pass"），但原位是 10.921 / 10.986 Å——分子形状对、位置完全错。这正是本次修复要消灭的误判。

### 4.3 平移探针（真实坐标上复现回归）

把上面运行产出的 `pose_1.sdf` 整体沿 x 平移 +10 Å（内部几何不变），再用**同一套代码**测量：

| 文件 | 原位 RMSD | 叠合 RMSD | D9 |
|---|---|---|---|
| `pose_1.sdf` 原样 | 0.454 Å | 0.380 Å | PASS |
| `pose_1_translated_10A.sdf` | **9.903 Å** | **0.380 Å（不变）** | **FAIL** |

修前用同一个文件会报 **0.380 Å / PASS**。

### 4.4 产物渲染核对（真数据跑真实消费点）

`/home/kui/qa-runs/d9-inplace/report.html` 的 `redocking` 段实际渲染（节选）：

```
<tr><td>1</td><td><b>0.454</b></td><td>✅ PASS</td><td class="muted">0.380</td></tr>
<tr><td>2</td><td><b>10.921</b></td><td>❌ FAIL</td><td class="muted">0.428</td></tr>
...
best-of-9（按 D9 原位口径取最优）：0.454 Å → PASS
参考配体：MK1:B。判据口径 inplace_calc_rms（原位 RMSD（不叠合，D9 判据），rdMolAlign.CalcRMS，不叠合）；
辅助口径 superposed_best_rms（叠合 RMSD（先最优叠合，仅辅助），rdMolAlign.GetBestRMS）
对刚体位移不敏感，不可用作 D9 结论。
```

audit Markdown 的 D9 段同样逐姿态列出并注明两个函数与"不可作为 D9 结论"。

---

## 5. 已知但**不在本任务范围**的问题：demo 路径产出非对映体

`/api/demo/example` 与 `server/tests/fixtures/indinavir_smiles.txt` 用的 indinavir SMILES 是**无立体描述**的：

```
CC(C)(C)NC(=O)C1CN(CCN1CC(CC(CC2=CC=CC=C2)C(=O)NC3C(CC4=CC=CC=C34)O)O)CC5=CN=CC=C5
```

它落地成 3D 时由 ETKDG 采样决定立体化学，实测得到的是 MK1 的**非对映体**：

| 分子 | InChIKey |
|---|---|
| 晶体 MK1:B（基准） | `CBVCZFGXHXORBI-PXQQMZJSSA-N` |
| demo SMILES（输入） | `CBVCZFGXHXORBI-UHFFFAOYSA-N`（立体未定） |
| demo SMILES 实际对接出的分子 | `CBVCZFGXHXORBI-RXFGFDIZSA-N`（非对映体，2 个手性中心与晶体相反） |

**本报告 §4 的实测刻意使用金标准立体 SMILES，因此 D9 数字是"真分子"的成绩。** demo 路径的这个缺陷不由本任务修复（属后续任务：demo SMILES 立体化学 / 配体身份校验）。在此之前，**任何人用 `/api/demo/example` 跑 D9，测的都是非对映体**——口径修好了，但分子还不是对的。

---

## 6. 验证结果

| 门禁 | 命令 | 结果 |
|---|---|---|
| crystal/rmsd 子集 | `PYTHONPATH=… pytest server/tests -q -k "crystal or rmsd"` | **21 passed, 1 skipped**（修前同命令 5 passed, 1 skipped） |
| **全量后端 pytest** | `PYTHONPATH=… DOCKSTUDIO_VINA=… pytest server/tests -q` | **553 passed, 25 skipped, 0 failed，957.63 s** ✅ |
| 收集用例总数 | `pytest --collect-only` | 562 → **578**（+16 条新用例：6 度量 + 10 接线） |
| lint | `ruff check server/dockstudio`（ruff 0.16.10） | **All checks passed!** |
| format | `ruff format --check server/dockstudio` | **53 files already formatted** |
| 真机 | §4 | 原位 0.454 Å PASS / 叠合 0.380 Å |

**零新增失败**：全量 0 failed（基线 `main` 上本机为 1 failed/497 passed 的既有 flaky 已在此轮环境不复现；本次运行本身 0 failed）。

---

## 7. 改动清单

分支 `fix/rmsd-inplace`（未合 main、未 push）：

| 文件 | 改动 |
|---|---|
| `server/dockstudio/analysis/crystal.py` | +215/−17：新增 `CrystalPoseRmsd`、`crystal_pose_rmsd`、`inplace_heavy_atom_rmsd`、`superposed_heavy_atom_rmsd`、`cocrystal_reference`、`pose_rmsd_series`、`best_pose_rmsd`、指标常量；`heavy_atom_rmsd` 保留并明确标注为叠合口径（deprecated as verdict） |
| `server/dockstudio/api/app.py` | `crystal_rmsd` 端点返回双口径；新增 `CrystalRmsdResponse` 响应模型（就近声明，避免动 `models.py`）；新增 `_job_cocrystal_rmsd_series` 供报告与审计共用；报告/审计端点传入 D9 数据 |
| `server/dockstudio/analysis/report.py` | +88：`redocking` 段（默认开启，可 `sections` 过滤）+ `_redocking_html` |
| `server/dockstudio/analysis/audit.py` | `report.md` D9 段 + `snapshot.json` 的 `cocrystal_rmsd` 块（`build_audit_markdown` 是 CLI/API 共享入口） |
| `server/tests/test_crystal_rmsd_metric.py` | 新增（6 条）：度量口径契约 |
| `server/tests/test_crystal_d9_wiring.py` | 新增（10 条）：API / 报告 / 审计 接线 |
| `docs/fixes/rmsd-inplace-measure_1hsg_d9.py` | 新增：真机测量脚本（本报告 §4 的产出方式） |
| `docs/fixes/rmsd-inplace-2026-11-18.md` | 本文件 |

`server/dockstudio/core/`、`prep/`、`apps/`、`resources/`、`scripts/`、`docs/qa/` **均未改动**。

> 范围说明：本任务的 inScope 是 `server/dockstudio/analysis`、`server/dockstudio/api`、`server/tests`、`docs/fixes`。因此 D9 响应模型就近声明在 `api/app.py`（该文件本就有 `CreateProjectRequest` / `PocketsRequest` / `FlexRequest` 等局部模型），而 D9 的 Markdown 段落在 `analysis/audit.py` 的 `build_audit_markdown` 里生成——CLI 的 `dockstudio audit` 本来就调用这个函数，所以**无需修改 `server/dockstudio/models.py` 或 `server/dockstudio/cli/`** 即可让 CLI 拿到 D9 段。

---

## 8. 概述本轮边界（明确没做的事）

| 事项 | 归属 |
|---|---|
| demo / fixture 的 indinavir SMILES 立体化学 | 后续任务（§5），本任务只记录事实 |
| 配体↔共晶配体身份校验 | 后续任务 |
| `DockParams.seed` 默认 `None` | 后续任务 |
| golden 集 RMSD 阈值收紧 | 后续任务 |
| 前端徽标改用原位值 | 前端不在 inScope；但 `rmsd` 字段语义已改，前端现有读取路径**自动**拿到正确的 D9 数字 |

## 9. 残余风险

1. **`heavy_atom_rmsd` 的语义在保存的旧数据里会"变味"**：历史 job 落盘的 `CrystalRmsd.rmsd` 是叠合值，新代码读出同一字段名会当成原位值。由于 `rmsd` 只由端点实时计算、不持久化在 job 结果里（`job.result.poses` 只有 affinity/rmsd_lb/rmsd_ub），实际不存在旧数据污染；若后续引入缓存，需带 `metric` 字段一并缓存。
2. **`symmetrizeConjugatedTerminalGroups`** 采用 RDKit 默认（`True`）。对羧酸/硝基等末端共轭基团会做对称校正，与文献常用实现一致；若将来要与某个特定对照工具逐位对齐，需要显式固定该参数。
3. 真机只测了 1HSG 一个靶点、一个种子（seed=42）。度量口径的正确性不依赖靶点（已由合成用例覆盖），但"该靶点该种子下 top1=0.454"这个具体数字不可外推到其他靶点/种子。