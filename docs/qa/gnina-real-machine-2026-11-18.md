# gnina 适配器真机验证报告（Tesla P4 · gnina 1.3.3）

- **日期**：2026-11-18
- **任务**：t2 / [t-gnina] —— 首次在真 GPU 上验证 DockStudio 的 gnina 适配器
- **性质**：只读验证。**未修改 `server/dockstudio/**`**，所有缺陷只记录不修复。
- **分支/工作树**：`qa/gnina-real` @ `/home/kui/wt/gpu`（基线 `main` = `c0f33cc`）
- **新增物**：本文件 + `server/tests/test_gnina_real_gpu.py`（opt-in 真机测试）

---

## 0. 结论摘要（TL;DR）

| 问题 | 结果 |
|---|---|
| GPU 路径是否真的能跑通？ | **能**。gnina 1.3.3 在 P4 上完成一次完整对接（28.5 s），nvidia-smi 全程可见 `engines/gnina` 进程，显存 8→202 MiB，GPU 利用率峰值 62% —— **没有静默回退 CPU**。 |
| CNN 打分能否解析？ | **能**，但来源只有 `out.pdbqt` 的 `REMARK CNNscore/CNNaffinity`。适配器为 stdout 表格写的那套解析在真机上是**死代码**，且一旦启用会**错列**（见 GN-02）。 |
| 原位重打分（`--score_only` / `--minimize`）可用吗？ | **不可用**。`dock()` 把任何配体一律改名成 `ligand.pdbqt`，导出的 `pose_N.sdf` 被当成 PDBQT 解析 → `Parse error ... Unknown or inappropriate tag`。生产里 Tier-2 的 "GNINA CNN rescore" 因此**永远失败且被 `except Exception: pass` 吞掉**（GN-08，本轮最严重）。 |
| 错误信息说人话吗？ | 一半一半。非法盒子/坏 PDBQT 有准确英文原因并可映射中文；**输入文件不存在时抛裸 `FileNotFoundError`**（GN-01），**未知打分函数时真实原因被截断**（GN-06），**空受体直接跑成功并返回 affinity=0.00 的"pose"**（GN-04）。 |
| `DOCKSTUDIO_GNINA_CPU=1` 回退可达吗？ | **可达**。命令追加 `--no_gpu`，无 GPU 时 `is_available()` 由"不可用"翻为"可用"，实测 CPU 运行 26.2 s 出正常 pose 与 CNN 值（但 `--no_gpu` 仍占 ~109 MiB 显存）。 |
| 文档承认"GPU 路径全是 mock 断言"是否属实？ | **属实**，且 mock 断言与真机有 6 处实质偏差（§4 表）。 |

---

## 1. 环境基线

| 项 | 值 |
|---|---|
| 主机 | `kui@100.98.223.80`（hostname `kui`，同时是 GitHub Actions 自托管 runner `ubuntu-p4`） |
| 内核 / CPU / 内存 | `Linux 7.0.0-38-generic` / 4 核 |
| GPU | Tesla P4，8192 MiB，驱动 580.178.04，CUDA 13.0（nvidia-smi 报告） |
| Python | `/home/kui/miniconda3/envs/ai-gpu/bin/python` = 3.11.15（rdkit 2026.03.6 / meeko 0.8.0 / torch 2.7.1+cu118） |
| gnina | `1.3.3 master:6fe1ce2  Built Jun 30 2026`（包装脚本 `resources/bin/linux-x86_64/gnina` → `exec /home/kui/engines/gnina`，`LD_LIBRARY_PATH` 指向 `~/engines/cudalib` + `~/engines/cudnn`） |
| 仓库基线 | `c0f33cc`（main），工作树 `/home/kui/wt/gpu` 分支 `qa/gnina-real` |
| 并发约束 | 所有真机运行均通过 `flock /home/kui/.dockstudio-heavy.lock` 串行化；每次运行前 `pgrep -af 'Runner.Worker\|pytest'` 记录并发状态 |
| 二进制定位 | `resources/bin/` 在 `.gitignore` 里（第 11 行），**新建工作树里没有 gnina 包装脚本**。本次已从主工作树复制 `gnina` + `vina` 到 `/home/kui/wt/gpu/resources/bin/linux-x86_64/`（该目录不入库）；未复制时 `is_available()` 会返回"gnina 二进制未找到：设置 DOCKSTUDIO_GNINA 或放入 resources/bin/linux-x86_64。" |

> 说明：多轮运行期间机器上还有队友的摸底 pytest / redock 批次在排队，本报告的每次 GPU 运行都在 `flock` 独占期内完成；`lock_*.log` 里保留了"取得锁时的并发进程快照"。

---

## 2. 测法与输入

**原则：不手敲 gnina 命令当结果。** 全部走 `dockstudio` 自己的代码路径：

```
server/tests/fixtures/1HSG.pdb
   └─ dockstudio.prep.receptor.prepare_receptor()        → receptor.pdbqt（链 A/B，0 warning）
server/tests/fixtures/indinavir_smiles.txt
   └─ dockstudio.prep.ligand.prepare_ligand_from_smiles() → ligand.pdbqt（45 重原子）
        └─ dockstudio.engines.gnina.GninaEngine().dock()  ← 被测对象
             └─ dockstudio.analysis.export.export_poses() → pose_1..N.sdf
```

盒子取自 1HSG 共晶配体 MK1 的质心（与 `test_e2e_1hsg.py` 同一套输入）：

```
center = (13.073, 22.467, 5.557)   size = (22.5, 22.5, 22.5)
```

探针脚本：`/home/kui/gnina_real_probe.py`、`/home/kui/gnina_probe2.py`（全文见附录 A）。

---

## 3. 证据

### 3.1 适配器生成的命令（逐字，来自 `dockstudio.engines.vina` 的 `running vina:` 日志）

```
running vina: /home/kui/wt/gpu/resources/bin/linux-x86_64/gnina --receptor receptor.pdbqt --ligand ligand.pdbqt --center_x 13.073 --center_y 22.467 --center_z 5.557 --size_x 22.5 --size_y 22.5 --size_z 22.5 --exhaustiveness 2 --num_modes 3 --cpu 3 --scoring vina --out out.pdbqt --seed 42
```

gnina 自己回显的 `Commandline:` 与此完全一致（除二进制绝对路径外）：

```
Commandline: /home/kui/engines/gnina --receptor receptor.pdbqt --ligand ligand.pdbqt --center_x 13.073 --center_y 22.467 --center_z 5.557 --size_x 22.5 --size_y 22.5 --size_z 22.5 --exhaustiveness 2 --num_modes 3 --cpu 3 --scoring vina --out out.pdbqt --seed 42
```

**命令拼装结论**：`--energy_range` 被成功剔除；`--cpu` / `--scoring vina` / `--seed` 都是 gnina 1.3.3 真实存在的开关（`gnina --help` 已逐条核对：`--cpu`、`--no_gpu`、`--score_only`、`--minimize`、`--scoring {ad4_scoring,default,dkoes_fast,dkoes_scoring,dkoes_scoring_old,vina,vinardo}` 全部存在）。**这一步 mock 断言是对的。**

### 3.2 真实 stdout（节选，`run/adapter_log.txt` 原样）

```
gnina v1.3.3 master:6fe1ce2   Built Jun 30 2026.
gnina is based on smina and AutoDock Vina.
Please cite appropriately.

Commandline: /home/kui/engines/gnina --receptor receptor.pdbqt --ligand ligand.pdbqt --center_x 13.073 --center_y 22.467 --center_z 5.557 --size_x 22.5 --size_y 22.5 --size_z 22.5 --exhaustiveness 2 --num_modes 3 --cpu 3 --scoring vina --out out.pdbqt --seed 42
Using random seed: 42

0%   10   20   30   40   50   60   70   80   90   100%
|----|----|----|----|----|----|----|----|----|----|
***************************************************
 | pose 0 | initial pose not within box

mode |  affinity  |  intramol  |    CNN     |   CNN
     | (kcal/mol) | (kcal/mol) | pose score | affinity
-----+------------+------------+------------+----------
    1      -11.19       -0.55       0.8610      8.779
    2      -10.47        1.34       0.7375      9.173
    3      -11.14       -0.26       0.4596      7.980
```

**注意两件事**：
1. 表头是**两行**、数据行是**五列**（`mode | affinity | intramol | CNN pose score | CNN affinity`）。适配器按四列写的正则见 GN-02。
2. 顺序不是能量序：`1 → -11.19`、`2 → -10.47`、`3 → -11.14`。gnina 默认 `--pose_sort_order=CNNscore`（`--help`：`How to sort docking results: CNNscore (default), CNNaffinity, Energy`），所以按 CNN score 0.8610 / 0.7375 / 0.4596 降序排列。见 GN-05。

另有一条 `| pose 0 | initial pose not within box`：配体初始构象不在盒内，属于 gnina 的**告警**而非错误，本次运行照常完成。

### 3.3 `out.pdbqt`（前 4 行真实内容）

```
MODEL 1
REMARK minimizedAffinity -11.194829
REMARK CNNscore 0.860957563
REMARK CNNaffinity 8.77859974
```

三个 REMARK 都在，且精度高于 stdout 表格。**适配器正是靠这三行拿到分值的**（`_AFF_RE` / `_CNN_SCORE_RE` / `_CNN_AFF_RE`）。

### 3.4 解析出的结构化值（`DockResult.poses`，`model_dump()` 原样）

| index | affinity (kcal/mol) | CNN score | CNN affinity (pKd) |
|---|---|---|---|
| 1 | **-11.194829** | 0.860958 | 8.778600 |
| 2 | **-10.4713345** | 0.737501 | 9.173482 |
| 3 | **-11.1351681** | 0.459579 | 7.979879 |

- `engine = "gnina"`，`engine_version = "1.3.3"`，`run_dir = /home/kui/gnina_probe_runs/dock/run`
- `rmsd_lb = rmsd_ub = None`（gnina 的 out.pdbqt 没有 vina 的 `REMARK VINA RESULT` 行，见 GN-12）
- 解析与独立复算一致：`parse_gnina_out(out.pdbqt, log)` 得到的 3 个 pose 与 `dock()` 返回值逐字段相同（该断言已固化为测试 `test_adapter_parse_matches_the_output_file`）。
- 导出的 SDF 可被 RDKit 读回，坐标零漂移：`pose_1.sdf` 首原子 `14.1730 24.0660 4.4220` == out.pdbqt 首原子 `14.173 24.066 4.422`。

### 3.5 "真上 GPU" 的证据（nvidia-smi 采样，200 ms 一次，共 118 个采样点 / 28.51 s）

```
distinct_compute_apps:
  873220, /home/kui/engines/gnina,   8 MiB
  873220, /home/kui/engines/gnina, 120 MiB
  873220, /home/kui/engines/gnina, 146 MiB
  873220, /home/kui/engines/gnina, 166 MiB
  873220, /home/kui/engines/gnina, 182 MiB
  873220, /home/kui/engines/gnina, 202 MiB
max_gpu_util_pct: 62
max_mem_used_mib: 205
```

即：gnina 进程确实持有 CUDA 上下文并在 P4 上跑了卷积网络打分，**不是**静默回退到 CPU。运行前后空闲态显存均为 0 MiB。

---

## 4. mock 断言 vs 真机差异清单

| # | 现有 mock/fixture 断言 | 真机事实 | 判定 |
|---|---|---|---|
| 1 | `REMARK minimizedAffinity` + `REMARK CNNscore` + `REMARK CNNaffinity`（`test_gpu_adapters.py::GNINA_OUT`） | 三行都存在，格式完全一致 | ✅ 吻合 |
| 2 | stdout 表格是**四列** `mode affinity CNNscore CNNaffinity`（`test_scoring_rigor.py::test_parse_gnina_stdout_mode_table`） | 1.3.3 是**五列**（多一列 `intramol`），且表头分两行 | ❌ 不吻合（GN-02） |
| 3 | 日志里会出现 `CNNscore 0.98 CNNaffinity -9.9` 这种"日志式"行（`GNINA_LOG`） | dock 模式**从不**打印该行；只有 `--score_only/--minimize` 的 stdout 会打印**冒号式** `CNNscore: 0.86097` / `CNNaffinity: 8.77868` | ❌ 半真（冒号式真实存在，空格式组合行不存在） |
| 4 | `cnn_affinity` 取负值（-9.9） | 真机是**正数 pKd**（7.98 / 8.78 / 9.17），与 `scoring_rigor.metric_directions["cnn_affinity"]=True`（越大越好）一致 | ❌ fixture 符号错（GN-10） |
| 5 | pose 列表按 affinity 升序（vina e2e 的通用契约） | gnina 按 CNN score 排序，affinity 非单调 | ❌ 不吻合（GN-05） |
| 6 | `dock()` 失败一律抛 `DockingError`（`base.py` 契约注释明写） | 输入文件缺失 → 裸 `FileNotFoundError`；空受体 → **不报错** | ❌ 不吻合（GN-01 / GN-04） |
| 7 | `--energy_range` 必须剔除；`--no_gpu` 存在；manifest 的 4 个 scoring 都被支持 | 全部成立 | ✅ 吻合 |
| 8 | `_binary()` 找不到 gnina 时给出中文提示 | 成立（但新工作树缺 `resources/bin` 时会命中，属于部署事实而非缺陷） | ✅ 吻合 |

**净结论**：`test_gpu_adapters.py::test_gnina_mock_binary_e2e` 这条"全流程"测试，跑的是一个会把 `GNINA_OUT` 原样写进 `--out` 的 shell 脚本 —— 它验证的是"适配器会不会读自己造的文件"，而不是"gnina 会不会产出这种文件"。真机一旦接入，第 2/4/5/6 条立刻不成立。

---

## 5. 缺陷清单（GN-01 … GN-12）

> 严重度：`blocker` 会直接产出错误科学结论；`high` 会让功能不可用或失败原因丢失；`medium` 影响可解释性/一致性；`low` 属于观感或未来风险。**全部只记录，不在本任务修复。**

### GN-08 [high] 原位重打分（`--score_only` / `--minimize`）在真机上完全不可用，且失败被静默吞掉

**现象**：把 `export_poses()` 导出的 `pose_1.sdf`（生产里 `screen_runner` Tier-2 CNN 复核用的正是它）喂给 `GninaEngine.dock(..., score_only=True)`：

```
DockingError: vina exited with code 1: ... | --- stderr --- | Parse error on line 2 in file "ligand.pdbqt": Unknown or inappropriate tag
PDBQT 文件解析失败：请重新准备受体/配体（meeko 输出可能不完整），并检查文件未被手工编辑。
```

**根因**：`VinaEngine.dock()` 第 339-344 行把任何配体一律 `shutil.copy2(ligand_pdbqt, work_dir / "ligand.pdbqt")` —— **扩展名被强行改成 `.pdbqt`**。gnina 按扩展名选解析器，于是 OpenBabel 把 SDF 当 PDBQT 读，第 2 行（含原子数的 counts 行）就炸。

**对照实验（证明不是 gnina 不会读 SDF，而是改名害的）**：同一条命令，仅把文件保持 `.sdf`：

```
$ gnina --receptor receptor.pdbqt --ligand ligand.sdf ... --score_only
Affinity: -11.19187 (kcal/mol)
CNNscore: 0.86097
CNNaffinity: 8.77867
→ exit 0, out.pdbqt 1 个 MODEL
```

**同样喂多模型 `out.pdbqt` 也失败**（`Parse error on line 92 ... Unknown or inappropriate tag`，第 92 行是 `ENDMDL`）；去掉 `ENDMDL` 的单模型 pdbqt 才能成功。

**生产影响（关键）**：`core/screen_runner.py:2291` 的 `top_pose_file` 来自 `_find_pose_file_in_dir()`，该函数优先返回 `<run_dir>/out.pdbqt`（多模型、带 `ENDMDL`）。于是 Tier-2 的：

```python
gnina_res = gnina_eng.dock(self.receptor_pdbqt, top_pose_file, ..., minimize=True, ...)
...
except Exception:
    pass          # ← screen_runner.py:2357-2358
```

**每一次都失败、每一次都被 `except Exception: pass` 吞掉**，候选表里 `cnn_score`/`cnn_affinity` 恒为空，报告不会出现任何"降级"标记。这正是团队第 2 轮要统一治理的"静默科学降级"的又一实例。

**最小复现**：

```python
from pathlib import Path
from dockstudio.analysis.export import export_poses
from dockstudio.engines.gnina import GninaEngine
from dockstudio.models import Box, DockParams

pose1 = export_poses(run_dir / "out.pdbqt", run_dir)[0]        # pose_1.sdf
GninaEngine().dock(receptor, pose1, box,
                   DockParams(exhaustiveness=1, num_modes=1, score_only=True),
                   run_dir / "score_only")
# → DockingError: Parse error on line 2 in file "ligand.pdbqt"
```

---

### GN-04 [high] 空受体被当成合法输入，返回 affinity=0.00 的"pose"

**现象**：把受体换成 0 字节文件，`dock()` **正常返回**：

```json
{"index": 1, "affinity": 0.0, "cnn_score": 0.269376576, "cnn_affinity": 4.09587908}
```

gnina 侧没有任何报错（网格里没有任何原子，能量恒等于 0）。适配器也不做任何受体完整性检查，于是"对接成功"这一结论被写进 `DockResult`。

**影响**：任何受体准备失败（空文件、只含 HETATM、meeko 静默产出 0 原子）都会变成"每个配体亲和力 0.0、CNN 值随机"的整批假结果，且没有任何降级标记。

**最小复现**：

```python
empty = tmp_path / "empty_receptor.pdbqt"; empty.write_text("")
res = GninaEngine().dock(empty, ligand, box, DockParams(exhaustiveness=1, num_modes=1), tmp_path / "run")
assert res.poses[0].affinity == 0.0     # 不抛异常
```

---

### GN-05 [high] pose 顺序是 CNN score 序，不是 affinity 序，但全仓把 `poses[0].affinity` 当 `top_affinity`

**事实**：gnina 默认 `--pose_sort_order=CNNscore`。实测 affinity 序列 `-11.19 / -10.47 / -11.14`，**非单调**；CNN score 序列 `0.8610 / 0.7375 / 0.4596` 严格递减。适配器从不重排。

**受影响调用点**（`grep 'poses\[0\]'`，均把第一个当成最优）：

```
api/app.py:364,626            core/queue.py:277,285        core/batch.py:143
core/project.py:30            core/screen_runner.py:1624,1767,1844,1998
cli/commands/dock.py:321      cli/commands/verify.py:531
```

后果：`ComparisonRow.top_affinity`、批量/筛选排名、`verify` 的 `reported_score` 比较，对 gnina 都是"CNN 最优位姿的 Vina 能量"，与 vina 的"能量最优位姿的 Vina 能量"**不同义**，跨引擎比较/共识排序被污染。

**可选修法（留给修复任务）**：给 gnina 追加 `--pose_sort_order 2`（Energy），或适配器返回前按 `affinity` 重排并在 `DockResult` 里保留原始顺序。

**最小复现**：`[p.affinity for p in res.poses] == sorted(...)` 在真机 3 个 pose 上就不成立。

---

### GN-01 [high] 输入文件缺失时抛裸 `OSError`，违反 `DockEngine` 契约，且没有中文提示

```
FileNotFoundError: [Errno 2] No such file or directory: '/home/kui/gnina_probe_runs/does_not_exist.pdbqt'
IsADirectoryError: [Errno 21] Is a directory: '/home/kui/gnina_probe_runs'      # 传目录也一样
```

异常 MRO：`FileNotFoundError → OSError → Exception`，**既不是 `DockingError` 也不是 `DockStudioError`**。而 `engines/base.py:55` 的契约写着"Must raise DockingError on failure, never leak raw subprocess errors"。API 路径的输入来自 registry 一般不会缺，但 CLI（`dockstudio dock --receptor ...`）与 `screen_runner` 直接传路径时会把英文 errno 甩给用户。

**最小复现**：`engine.dock(receptor, Path("does_not_exist.pdbqt"), box, params, work)`。

---

### GN-02 [medium] stdout 表格解析：五列被当四列（一旦启用就会错列），且表头检测在 1.3.3 上永远失败

两件事叠加：

1. **表头正则失效**：`_GNINA_TABLE_HEADER_RE = r"CNN\s*(?:pose\s*)?score|CNN\s*affinity"` 要求 `CNN` 与 `pose score` 之间只有空白。1.3.3 的表头是**两行**：

   ```
   mode |  affinity  |  intramol  |    CNN     |   CNN
        | (kcal/mol) | (kcal/mol) | pose score | affinity
   ```

   中间隔着 `| (kcal/mol) | (kcal/mol) |`，`\s*` 跨不过去 → `has_gnina_table_header == False`（实测）。所以四列分支在真机上是**死代码**，真机走的是 vina 兼容分支，只拿到 affinity，CNN 两列被**静默丢弃**。

2. **就算表头命中也会错列**：数据行是五列，四列正则 `(\d+)\s+(float)\s+(float)\s+(float)` 取前四个数字。用真实数据行验证：

   ```
   输入行: "    1      -11.19       -0.55       0.8610      8.779"
   groups: ('1', '-11.19', '-0.55', '0.8610')
   → cnn_score = -0.55   （其实是 intramol 列）
   → cnn_affinity = 0.8610（其实是 CNN pose score）
   ```

   单行表头（老版本 gnina 写法）下的完整解析结果确认了这一点：`cnn_score=-0.55 / cnn_affinity=0.861`，全是**错列**。

**最小复现**：

```python
from dockstudio.engines.gnina import parse_gnina_out
log = Path("adapter_log.txt").read_text()          # 真机 stdout
parse_gnina_out(Path("/tmp/empty.pdbqt"), log)     # → 3 个 pose，cnn_* 全为 None
```

---

### GN-03 [medium] `out.pdbqt` 缺失时静默降级：CNN 两个字段变 `None`，无警告无标记

`parse_gnina_out` 在 `out.pdbqt` 不存在/不可解析时回退到 stdout 表格，返回的 `Pose.cnn_score/cnn_affinity` 为 `None`，**没有任何 `log.warning`、没有 Pose 级标记、没有写进报告**。结合 GN-02（表头检测失效），这条回退路径在 1.3.3 上必然产出"有 affinity、无 CNN"的 pose。下游 `screen_runner.py:2355` 直接把 `None` 写进候选行。

### GN-06 [medium] `DockingError` 截断丢掉真实原因（`_key_lines(..., 5)` 的副作用）

未知打分函数时：

```
DockingError: vina exited with code 1: dkoes_fast | dkoes_scoring | dkoes_scoring_old | vina | vinardo
```

真正的报错（`ERROR: could not find scoring function: not_a_function`）在 stdout 前部，被"取最后 5 行"截掉了，用户只看到一串合法函数名。同理，很多 gnina 失败会把 `Please cite appropriately.` / `Commandline: ...` 当成"关键信息"顶到最前面（见 GN-08 的报文）。

### GN-07 [low] 盒子没有前置校验：0/负数直接送进 gnina，60 Å 也照跑

- `Box(center, size)` 无任何 validator；`size=(0,0,0)` / `(-5,-5,-5)` 会一路传到子进程，靠 gnina 报 `Usage error: Search space dimensions should be positive`（**这一条被正确转成 `DockingError`，且没有被误判成"盒子超限"** —— 见 §6，是好消息）。
- `size=60³` 不报错、正常出结果（`-11.18`），而 `vina._ERROR_KB` 的中文提示写着"建议各边 ≤ 30 Å" —— 该建议在 gnina 侧从未被执行，也没有任何 warning。

### GN-09 [low] `--score_only` 输出里出现两个 `MODEL 1`，导致 `Pose.index` 重复

真机 `--score_only` 的 `out.pdbqt` 结构：

```
1: MODEL 1
2: REMARK minimizedAffinity -11.1918707
3: REMARK CNNscore 0.860967219
4: REMARK CNNaffinity 8.77867508
5: MODEL 1                     ← 又一个 MODEL 1
6: REMARK minimizedAffinity -11.194829
...
96: ENDMDL                      ← 只有一个 ENDMDL
```

`parse_gnina_out` 于是返回两个 `index=1` 的 pose（`-11.1918707` 和 `-11.194829`）。前端/导出若用 `index` 当 key 会重复。

### GN-10 [low] mock fixture 与真机格式偏离（详见 §4 表第 2/3/4 条）

`test_scoring_rigor.py` 里三条 `test_parse_gnina_*` 与 `test_gpu_adapters.py::GNINA_OUT/GNINA_LOG` 构造了真机不会产出的文本（四列表、空格式组合行、负的 CNN affinity）。它们保住了"正则能匹配字符串"，但没有保住"能匹配 gnina"。

### GN-11 [low] 运行日志与报错文案把引擎名写死成 `vina`

`vina.py` 里 `log.info("running vina: ...")`、`DockingError("vina exited with code ...")`、`"vina timed out after ..."`、`"vina produced no output file"` 对 gnina / smina / Uni-Dock 全部复用。用户看到"vina exited with code 1"时其实跑的是 gnina，排查会被带偏。

### GN-12 [low] gnina 的 `rmsd_lb/rmsd_ub` 永远为 `None`

gnina 的 out.pdbqt 没有 `REMARK VINA RESULT`（那是 vina 的行），解析器也不会从别处补。于是 UI/报告里"RMSD l.b./u.b."两列对 gnina 恒空，而 vina 的通用契约是"poses 按 affinity 排序 + 带 RMSD"。若后续要复用 redock RMSD 判据，必须自己算（不能读这两个字段）。

---

## 6. 失败路径实测矩阵（"错误信息说人话吗、会不会误判成盒子超限"）

| 场景 | 适配器抛出 | 中文提示 | 判定 |
|---|---|---|---|
| 配体文件不存在 | `FileNotFoundError`（裸 OSError） | 无 | ❌ GN-01 |
| 受体文件不存在 | `FileNotFoundError`（裸 OSError） | 无 | ❌ GN-01 |
| 配体路径是目录 | `IsADirectoryError`（裸 OSError） | 无 | ❌ GN-01 |
| 配体是垃圾文本 | `DockingError` | ✅"PDBQT 文件解析失败：请重新准备受体/配体…" | ✅ 分类正确 |
| 受体是 0 字节 | **不抛异常**，返回 `affinity=0.00` | 无 | ❌ GN-04 |
| 盒子 size=(0,0,0) | `DockingError`，含 `Usage error: Search space dimensions should be positive` | 无中文 | ✅ 没误判成"超限"，但缺中文（GN-07） |
| 盒子 size=(-5,-5,-5) | 同上 | 无中文 | 同上 |
| 未知打分函数 | `DockingError: vina exited with code 1: dkoes_fast \| dkoes_scoring \| ...` | 无 | ❌ 真实原因被截断（GN-06） |
| 盒子 60³ | **不抛异常**，正常出结果 | 无 warning | ⚠️ GN-07（"≤30 Å"建议未执行） |
| 配体初始位置在盒外（6³ 盒） | **不抛异常**，返回 `affinity=-0.21` | 无 | ⚠️ 静默产出无意义结果 |

**关于"误判成盒子超限"**：本轮 10 个失败场景里**没有一例**被误标成"搜索盒设置超限"。`vina._ERROR_KB` 里那条 `box.*exceed|too large|search space is too|...` 的规则对 gnina 的报文基本不会被触发（gnina 用 `Usage error: Search space dimensions should be positive`，不匹配）。风险点反过来：**该报"盒子有问题"的场合（6³ 盒导致无意义结果）反而什么都不说**。

---

## 7. `DOCKSTUDIO_GNINA_CPU=1` 回退路径

**(a) 命令层面**：环境变量设置后 `_build_cmd` 追加 `--no_gpu`，逐字命令：

```
running vina: /home/kui/wt/gpu/resources/bin/linux-x86_64/gnina --receptor receptor.pdbqt --ligand ligand.pdbqt --center_x 13.073 --center_y 22.467 --center_z 5.557 --size_x 22.5 --size_y 22.5 --size_z 22.5 --exhaustiveness 1 --num_modes 1 --cpu 3 --scoring vina --out out.pdbqt --no_gpu --seed 42
```

**(b) 可用性门控**（无需 GPU 的纯逻辑验证）：

| 状态 | `is_available()` |
|---|---|
| GPU 在位，env 未设 | `(True, None)` |
| GPU 在位，env=1 | `(True, None)` |
| GPU 探测返回 False，**env=1** | `(True, None)` ← 逃生门有效 |
| GPU 探测返回 False，env 未设 | `(False, '…no NVIDIA GPU；设置 DOCKSTUDIO_GNINA_CPU=1 可强制 CPU（--no_gpu，极慢）')` |

**(c) 真实 CPU 运行**：`exhaustiveness=1, num_modes=1` 用时 **26.18 s** 完成，返回 `affinity=-11.2035875 / CNN score 0.8004 / CNNaffinity 8.575`，out.pdbqt 三个 REMARK 齐全 —— 回退路径**确实可达且能出正确结果**。

**(d) 一个诚实的观察**：`--no_gpu` 运行期间 nvidia-smi 仍报告该 gnina 进程占用 **~109 MiB** 显存、GPU 利用率峰值仅 2%。也就是说这个构建即使 `--no_gpu` 也会建 CUDA 上下文（打分在 CPU 上做，但"完全不碰 GPU"并不成立）。若未来要在无 CUDA 库的机器上跑，需要另行验证。

---

## 8. 建议的修复顺序（留给后续任务，本轮不改代码）

1. **GN-08**：`dock()` 保留输入配体的真实格式（按内容 sniff 或按扩展名保留文件名），或统一"预处理成 PDBQT 再喂 gnina"；同时把 `screen_runner.py` 的 `except Exception: pass` 换成"写降级台账 + 候选行标记"。
2. **GN-04**：`dock()` 前置校验受体 PDBQT 至少含 1 个 `ATOM`，否则 `DockingError`；对 `affinity == 0.0` 增加合理性告警。
3. **GN-05**：给 gnina 追加 `--pose_sort_order 2`（Energy），或返回前按 affinity 重排。
4. **GN-01**：`dock()` 入口校验 `is_file()`，缺失即 `DockingError("配体文件不存在：…")`。
5. **GN-02/GN-03**：按 1.3.x 的五列 + 两行表头重写 stdout 解析；回退时打 `HIGH_*_FALLBACK` 标记。
6. 其余（GN-06/07/09/10/11/12）按观感/一致性排期。

---

## 附录 A：探针脚本

- `/home/kui/gnina_real_probe.py` —— 主流程（prep → dock → 失败矩阵 → CPU 回退），输出 `evidence.json`
- `/home/kui/gnina_probe2.py` —— 原位重打分（SDF/PDBQT 两种输入）、错误文案、CPU 回退，输出 `evidence2.json`
- 现场产物：`/home/kui/gnina_probe_runs/`（`dock/run/{out.pdbqt,adapter_log.txt,pose_*.sdf}`、`fail_*/`、`sdf_*/`、`cpu_fallback/`、`lock_*.log`）

## 附录 B：opt-in 真机测试

`server/tests/test_gnina_real_gpu.py`（12 个用例，默认全部 skip）：

```
$ pytest tests/test_gnina_real_gpu.py -q
12 skipped

$ DOCKSTUDIO_REAL_GPU=1 DOCKSTUDIO_GNINA=…/gnina pytest tests/test_gnina_real_gpu.py -q -rxX
（结果见 §9）
```

设计要点：
- 默认 `skip`，只有 `DOCKSTUDIO_REAL_GPU=1` **且** gnina 可解析时才运行；CI 无卡环境零影响。
- 一条 `dock_run`（module 级 fixture）跑一次真机对接（≈30 s），供多数断言复用。
- GN-01/GN-04/GN-05/GN-08 的用例标成 `xfail(strict=False)`：现在 xfail，修复后变 XPASS 而不会让套件变红 —— 既是"缺陷快照"，又是"修复验收尺子"。
- GN-02（五列错列）用 §5 的最小复现留在报告里而**不入测试**，因为它需要一段来自真机的 stdout 才能复现，写进测试会变成"用 fixture 断言 fixture"的老毛病。

---

## 9. 真机测试运行结果

**默认（CI / 无卡机器）—— 必须全 skip：**

```
$ cd /home/kui/wt/gpu/server && python -m pytest tests/test_gnina_real_gpu.py -q
12 skipped, 1 warning in 0.04s
```

**真机（Tesla P4，`DOCKSTUDIO_REAL_GPU=1`）—— `real_gpu_test.log` 原样：**

```
$ DOCKSTUDIO_REAL_GPU=1 DOCKSTUDIO_GNINA=/home/kui/wt/gpu/resources/bin/linux-x86_64/gnina \
      python -m pytest tests/test_gnina_real_gpu.py -q -rxX
.....xxx.x..                                                             [100%]

slowest durations:
28.43s call     ...::test_empty_receptor_is_rejected
27.73s setup    ...::test_real_gpu_dock_returns_sane_poses
 0.78s call     ...::test_in_place_rescore_accepts_an_exported_pose
 0.53s call     ...::test_zero_size_box_raises_docking_error
 0.03s call     ...::test_exported_pose_keeps_coordinates

XFAIL ...::test_poses_are_sorted_best_first_by_affinity - GN-05: gnina sorts poses by CNN score (--pose_sort_order=CNNscore), the adapter never re-sorts by affinity; every call site uses poses[0].affinity
XFAIL ...::test_missing_ligand_raises_docking_error - GN-01: shutil.copy2's FileNotFoundError escapes dock() instead of DockingError
XFAIL ...::test_empty_receptor_is_rejected - GN-04: gnina happily docks against an empty receptor and reports affinity 0.00
XFAIL ...::test_in_place_rescore_accepts_an_exported_pose - GN-08: dock() renames every ligand to ligand.pdbqt, so an exported pose_1.sdf is parsed as PDBQT and gnina fails with 'Unknown or inappropriate tag'

8 passed, 4 xfailed, 1 warning in 57.67s
```

- **8 passed** 覆盖：真机出 pose（affinity 量级/符号合理）、`--no_gpu` 不出现（无静默 CPU 回退）、CNN 两个字段从 `out.pdbqt` 解析出且 `cnn_affinity > 0`、适配器解析与独立复算逐字段一致、导出 SDF 坐标零漂移、零尺寸盒抛 `DockingError` 且不含"搜索盒设置超限"、CPU 逃生门的命令层面与可用性层面。
- **4 xfailed** 就是 GN-01 / GN-04 / GN-05 / GN-08 的"缺陷快照"，全部 `strict=False`：现状 xfail，修复后 XPASS 而不会把套件弄红。
- 新文件本身通过 `ruff check` 与 `ruff format --check`（ruff 0.16.10，`line-length=100`）。
