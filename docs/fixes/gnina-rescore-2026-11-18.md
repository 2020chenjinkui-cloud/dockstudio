# GNINA 原位重打分修复报告（GN-08 / GN-03 / GN-04）

- **日期**：2026-11-18
- **任务**：t8 —— 修掉 Tier-2 GNINA 原位重打分完全不可用的缺陷，并消灭该路径上的静默降级
- **分支**：`fix/gnina-rescore` @ `/home/kui/wt/gnina-fix`（基线 `main` = `852d26c`；当前 tip `9397ecf`，9 文件 +1183 / −49）
- **未合 main、未 push**（按约定由 integrator 串行处理）
- **前置证据**：`docs/qa/gnina-real-machine-2026-11-18.md`（t2 的真机取证，GN-01…GN-12）

---

## 0. 摘要

| 缺陷 | 修前（真机/回归用例实测） | 修后（真机实测） |
|---|---|---|
| **GN-08** 原位重打分不可用 | `pose_1.sdf` 被改名 `ligand.pdbqt` → `Parse error on line 2 ... Unknown or inappropriate tag`；多模型 `out.pdbqt` → `Parse error on line 92`（ENDMDL） | `--ligand ligand.sdf`（SDF 保留扩展名）；多模型 → 只留 MODEL 1、去 ENDMDL，`--ligand ligand.pdbqt` 正常打分。P4 实测 `cnn_score=0.828008473 / cnn_affinity=8.72152042`（SDF 输入）、`0.827973425 / 8.72327232`（out.pdbqt 输入） |
| GN-08 静默吞异常 | `screen_runner.py` 里 `except Exception: pass`，候选行 `alerts=[]` | 候选行 `alerts` 出现 `GNINA_RESCORE_UNAVAILABLE:<真实原因>`；生产路径 `_run_tier3` 端到端仍为 `alerts=[]`（成功时），失败时留痕（负向对照 D） |
| **GN-03** CNN 缺失静默变 None | 无 warning、无标记 | `parse_gnina_out` 打 warning；`DockResult.log` 追加 `DEGRADED GNINA_CNN_MISSING:…`；Tier-3 候选行同样打 `GNINA_CNN_MISSING:` |
| **GN-04** 空受体伪成功 | 0 字节受体 → `affinity=0.0` + `cnn_score=0.2694` 的"pose"，视为成功 | 校验层拒绝：`DockingError: 受体文件不含任何原子记录（ATOM/HETATM）：… 空受体无法用于对接，重新准备受体后再试。`（真机回归用例由 fail→pass） |

顺带的收益：真机测试套件从 57.2 s 降到 **32.1 s**——空受体用例不再花 28 s 对着空气"对接"。

---

## 1. 缺陷与修复对照

| # | 现象（t2 真机取证） | 根因 | 修复落点 |
|---|---|---|---|
| GN-08a | SDF 位姿喂给 `--score_only/--minimize` 必炸 | `vina.py:339-344` 把任何配体 `shutil.copy2` 成 `ligand.pdbqt`，引擎按扩展名选解析器 | `VinaEngine._stage_ligand()`：扩展名保持（`.sdf/.mol2/.pdb/.mol/.pdbqt`），未知后缀回落 `.pdbqt` |
| GN-08b | 多模型 `out.pdbqt` 在 ENDMDL 处炸 | gnina PDBQT 解析器不接受 `ENDMDL`，只吃单 model | `GninaEngine._prepare_staged_ligand()`：MODEL 1 + 去 ENDMDL（取舍见 §6） |
| GN-08c | 失败被 `except Exception: pass` 吞掉，CNN 列恒空 | `screen_runner.py:2457-2458` | 抽出 `screen_runner.gnina_cnn_rescore()`，**任何**不成功都返回结构化标记；`_run_tier3` 把标记写进候选行 `alerts` |
| GN-03 | CNN 字段缺失静默 None | 解析器与适配器都不检查 | `parse_gnina_out()` warning + `GninaEngine._degradation_notes()` → `DockResult.log` 追加 `DEGRADED GNINA_CNN_MISSING:…` |
| GN-04 | 空受体返回 `affinity=0.0` 伪 pose | 无任何输入校验 | `engines/base.py::validate_receptor_file()`（带 `lru_cache` 的原子计数），`VinaEngine.dock()` 启动引擎前调用（smina/gnina 走该入口）；**unidock 的 `dock_batch()` 是独立入口，已单独接线**（T9-F2）。`BaseEngine.dock`/`dock_batch` 是协议/兜底，不会自动调用——不要写成「四个引擎全部继承」 |

---

## 2. 代码改动（9 文件）

| 文件 | 改动 | 为什么 |
|---|---|---|
| `server/dockstudio/engines/base.py` | +105：`GNINA_RESCORE_UNAVAILABLE` / `GNINA_CNN_MISSING` / `marker()` / `log_line()`（降级标记）；`SUPPORTED_LIGAND_SUFFIXES`；`validate_receptor_file()`（ATOM/HETATM 计数，按 `(path, mtime, size)` 记忆化） | 标记必须是**纯字符串**才能不改任何持久化 schema 就穿过 `alerts → CSV → CLI → checkpoint.jsonl → DockResult.log`（旧 checkpoint 仍可读）；输入卫生放协议层，所有适配器共享，HTVS 里每个配体都会调用必须便宜 |
| `server/dockstudio/engines/vina.py` | +46：新增 `_stage_ligand()` / `_prepare_staged_ligand()` / `_degradation_notes()` 三个钩子；`dock()` 改为"校验受体 → 按真实格式暂存配体 → 解析后收集降级说明" | 把"配体到底是什么格式"变成显式判断（扩展名 + 引擎钩子），而不是一律改名 |
| `server/dockstudio/engines/gnina.py` | +69：`_prepare_staged_ligand()`（MODEL 1 / 去 ENDMDL）、`_degradation_notes()`（CNN 缺失）、`parse_gnina_out()` warning、模块 logger | gnina 的两个真实解析怪癖只在 gnina 里处理，不污染 vina |
| `server/dockstudio/core/screen_runner.py` | +121：新增模块级 `gnina_cnn_rescore()`；`_run_tier3` 改为调用它并把标记塞进 `alerts`；GNINA 可用性探测由 `pass` 改为记日志 | 把原来内联在 `_run_tier3` 里的 `try/except pass` 抽成**可单测**的函数，标记语义集中在一处 |
| `server/tests/test_gnina_rescore_fix.py` | **新增 494 行 / 15 用例**（不依赖 GPU） | 先红后绿（§3） |
| `server/tests/test_gnina_real_gpu.py` | 2 个 `xfail` 转正式断言（原位重打分、空受体），并加 `cnn_affinity`/无降级标记断言 | GN-04/GN-08 已修，快照变验收 |
| `server/tests/test_gpu_adapters.py` | mock e2e 的受体从 `"x"` 换成含 ATOM 记录的 PDBQT | 否则新校验会让该用例**静默 skip**（实测发生），等于偷偷丢覆盖 |
| `server/tests/conftest.py` | +12：把本工作树的 `server/` 插到 `sys.path` 首位 | 见 §7.3：editable 安装的 meta-path finder 指向 `/home/kui/dockstudio`，从仓库根跑 pytest 会**测到 main 而不是被测分支**（t5-F1 的同一个坑） |

---

## 3. 先红后绿

失败用例先写在**未修**的 `fix/gnina-rescore`（= `852d26c`）上跑，4 条全红，且失败原因就是生产现场那两句：

```
$ cd /home/kui/wt/gnina-fix/server && python -m pytest tests/test_gnina_rescore_red_probe.py -q
E  dockstudio.models.DockingError: vina exited with code 1: --- stderr --- | Parse error on line 2 in file "ligand.pdbqt": Unknown or inappropriate tag
E  dockstudio.models.DockingError: vina exited with code 1: --- stderr --- | Parse error on line 92 in file "ligand.pdbqt": Unknown or inappropriate tag
E  Failed: DID NOT RAISE DockingError          # 空受体伪成功
E  AssertionError: []                          # 候选行 alerts 里没有降级标记（异常被吞）
4 failed, 1 warning in 5.42s
```

四个用例分别对应：SDF 位姿重打分、多模型 `out.pdbqt` 重打分、空受体、`_run_tier3` 标记。

同一份探针在修完后 **4 passed**；正式测试文件 `tests/test_gnina_rescore_fix.py` **15 passed**（1.81 s）。探针是临时文件，验证完即删除，留下的正式用例覆盖同样断言：

```
$ pytest server/tests/test_gnina_rescore_fix.py -q     →  15 passed, 1 warning in 1.81s
```

其中 3 条是"复现 + 反证"组合，值得单独点出：

- `test_in_place_rescore_accepts_an_exported_sdf_pose`：不仅要求成功，还从子进程 `argv.txt` 断言 `ligand.sdf` 在、`ligand.pdbqt` 不在 —— 直接锁死"改名"这件事。
- `test_multimodel_out_pdbqt_is_reduced_to_one_model`：断言暂存文件里 `MODEL` 只有 1 个、没有 `ENDMDL`、且留下的是 MODEL 1 的 `minimizedAffinity`。
- `test_tier3_candidate_carries_the_rescore_degradation_marker`：走真实 `ScreeningRunner._run_tier3`，断言 `alerts` 出现 `GNINA_RESCORE_UNAVAILABLE:`，并进一步断言它出现在 `screening_results.csv` 的 `Alerts` 列。

---

## 4. 真机证据（Tesla P4 8GB / gnina 1.3.3 / driver 580.178.04 / CUDA 13.0）

跑法（`flock` 独占，`DOCKSTUDIO_GNINA` 指向主工作树的包装脚本）：

```
$ flock /home/kui/.dockstudio-heavy.lock -c \
    "DOCKSTUDIO_GNINA=/home/kui/dockstudio/resources/bin/linux-x86_64/gnina \
     python /home/kui/gnina_rescore_gpu_probe.py"
[env] gnina available=True reason=None version=1.3.3
[prep] reusing /home/kui/gnina_probe_runs/dock/run/out.pdbqt
[prep] exported poses: ['pose_1.sdf', 'pose_2.sdf', 'pose_3.sdf']
```

探针调用的是 **DockStudio 自己的代码**（`screen_runner.gnina_cnn_rescore` / `ScreeningRunner._run_tier3` / `analysis.export.export_poses`），没有手敲 gnina 命令当结果。

### 4.1 A 组：SDF 位姿（修前必炸的那条输入）

命令原样（适配器打印的 `running vina:` 日志）：

```
running vina: /home/kui/dockstudio/resources/bin/linux-x86_64/gnina --receptor receptor.pdbqt --ligand ligand.sdf --center_x 13.073 --center_y 22.467 --center_z 5.557 --size_x 22.5 --size_y 22.5 --size_z 22.5 --exhaustiveness 1 --num_modes 1 --cpu 4 --scoring vina --out out.pdbqt --minimize
```

注意 `--ligand ligand.sdf`：**SDF 扩展名一路保持到子进程**。暂存目录里只有 `ligand.sdf`。

真实输出（`out.pdbqt` 前 5 行）：

```
MODEL        2
REMARK  Name =
REMARK minimizedAffinity -11.1307526
REMARK minimizedRMSD 0.352880895
REMARK CNNscore 0.828008473
```

结构化结果：`cnn_score = 0.828008473`，`cnn_affinity = 8.72152042`，`degradation_marker = None`。

GPU 证据：4.73 s，`nvidia-smi --query-compute-apps` 抓到 `923644, /home/kui/engines/gnina, 112→194 MiB`，显存峰值 **225 MiB**。

### 4.2 B 组：多模型 `out.pdbqt`（生产的真实输入形态）

`_find_pose_file_in_dir()` 优先返回 `<run_dir>/out.pdbqt`（3 个 MODEL），所以生产上传进去的就是这个文件。

```
running vina: …/gnina --receptor receptor.pdbqt --ligand ligand.pdbqt … --minimize
staged ligand: models=1  endmdl=False   (源文件 models=3)
→ cnn_score = 0.827973425, cnn_affinity = 8.72327232, marker = None
GPU: 1.28 s, 峰值 219 MiB, compute-app pid 923740
```

A 与 B 的 CNN 值有 4e-5 级差异，来源是 SDF 只有 4 位坐标精度而 PDBQT 是全精度 —— 属预期，不影响结论。

### 4.3 C 组：Tier-3 生产路径端到端

`ScreeningRunner(config.enable_gnina_rescore=True).__tier3__`，`best_pose_path` 就是上面那个 `out.pdbqt`：

```
[C] candidate cnn_score=0.827973425 cnn_affinity=8.72327232 alerts=[]
GPU: 13.48 s, 峰值 225 MiB, compute-app pid 923832
CSV: … CNN_Score=0.8280, CNN_Affinity=8.7233 …
```

也就是说：**候选表的 CNN 两列在真 GPU 上第一次拿到了真实数值**（修前恒为空）。

### 4.4 D 组：负向对照（真失败时必须留痕）

给一个垃圾 SDF 当位姿：

```
[D] broken ligand -> marker=GNINA_RESCORE_UNAVAILABLE:gnina output contained no binding poses:
gnina v1.3.3 master:6fe1ce2 Built Jun 30 2026. | gnina is based on smina and AutoDock Vina. |
Please cite appropriately. | Commandline: /home/kui/engines/gnina --receptor receptor.pdbqt
--ligand ligand.sdf …
```

标记里带的是**真实引擎文本**，不是我们编的字符串；同一路径下 `cnn_score/cnn_affinity` 为 `None`（不再伪装成数值）。

### 4.5 真机 pytest

```
$ cd /home/kui/wt/gnina-fix && DOCKSTUDIO_REAL_GPU=1 DOCKSTUDIO_GNINA=…/gnina \
    python -m pytest server/tests/test_gnina_real_gpu.py -q -rxX
10 passed, 2 xfailed, 1 warning in 32.07s
XFAIL test_poses_are_sorted_best_first_by_affinity  - GN-05（本轮不修）
XFAIL test_missing_ligand_raises_docking_error      - GN-01（本轮不修）
```

修前同一套件是 `8 passed, 4 xfailed`（GN-04/GN-08 还是 xfail）；现在这两条变成正式断言并通过。

---

## 5. 降级标记的四个可验证落点

标记格式：`GNINA_RESCORE_UNAVAILABLE:<原因>` / `GNINA_CNN_MISSING:<原因>`（`engines/base.py::marker`，单行、截断到 240 字符，CSV 单元格安全）。

| 落点 | 代码位置 | 怎么验证 |
|---|---|---|
| 候选行 `alerts` | `screen_runner._run_tier3` → `ScreeningCandidateRow(alerts=[…])` | `test_tier3_candidate_carries_the_rescore_degradation_marker` / `..._cnn_missing_marker` |
| `screening_results.csv` 的 `Alerts` 列 | `export_screening_results_csv()`（`;` 连接） | 同上用例断言 `GNINA_RESCORE_UNAVAILABLE in csv_text`；CSV 表头含 `Alerts` |
| `checkpoint.jsonl`（审计/续跑） | `CheckpointManager.record_tier3(lig.id, cand.model_dump())` | `test_candidate_row_keeps_alert_roundtrip`（`model_dump()` 往返不丢） |
| `DockResult.log`（单分子对接 / CLI / API 响应） | `VinaEngine.dock()` 末尾 `DEGRADED <标记>` 追加 | `test_cnn_missing_result_carries_a_degradation_marker` |

另外 `cli/commands/screen.py:422-426` 本来就把 `alerts` 渲染成候选表徽标，所以终端也能看到。

**没有**把"gnina 没装"这一部署事实盖到每一行上：`_run_tier3` 的可用性探测现在会 `log.warning`，但只有**真正尝试过**的重打分才在候选行留标记。理由：没装 gnina 是用户可从 `/api/engines` 直接看到配置事实，而"装了、跑了、失败了却没人知道"才是科学降级。

---

## 6. 格式处理与取舍

1. **SDF / MOL2 / PDB / MOL 输入**：暂存名保持原扩展名（`ligand.sdf` 等），绝不改写成 `.pdbqt`。
2. **多模型 PDBQT**（`out.pdbqt`）：只保留 MODEL 1，并删除 `ENDMDL`。取舍理由：
   - `--score_only` / `--minimize` 的语义就是"给一个位姿、就地打分/优化"，多模型本身没有意义；
   - Tier-2 传进来的 `out.pdbqt` 的第一个 model 正是本次要复核的最优位姿；
   - 若改成"逐 MODEL 各跑一次"，会让一次 Tier-3 复核变成 N 次 GPU 调用（N=num_modes，默认 9），成本 ×9 而科学价值不变。
   - 该行为只在 `.pdbqt` 且文件里出现 `MODEL`/`ENDMDL` 时触发；meeko 产出的普通单模型配体 byte-for-byte 不变（`test_plain_single_pose_pdbqt_is_left_untouched`）。
3. **未知后缀**（如 `.txt`）：仍按 `.pdbqt` 暂存，让引擎报解析错并由 `_ERROR_KB` 翻译成中文提示（保持既有行为）。

---

## 7. 门禁与不回归

### 7.1 三条 verify 命令

```
$ cd /home/kui/wt/gnina-fix && python -m pytest server/tests/test_gpu_adapters.py server/tests/test_scoring_rigor.py -q
35 passed, 2 skipped, 1 warning in 1.33s

$ cd /home/kui/wt/gnina-fix && DOCKSTUDIO_REAL_GPU=1 \
    python -m pytest server/tests/test_gnina_real_gpu.py -q
10 passed, 2 xfailed, 1 warning in 32.21s
# 不设 DOCKSTUDIO_REAL_GPU 时该文件整体 skip（12 skipped）。
# gnina 二进制从 resources/bin/linux-x86_64/gnina 解析（该目录在 .gitignore 里，
# 本工作树已从主工作树复制；详见 §7.4）。

$ cd /home/kui/wt/gnina-fix && ruff check server/dockstudio && ruff format --check server/dockstudio
All checks passed!
53 files already formatted
```

（`ruff` = `/home/kui/miniconda3/envs/ai-gpu/bin/ruff`，0.16.10，沿用 t4 锁定的显式 `select`。）

### 7.2 全量 pytest（引擎可用条件，flock 内独占）

全量共跑四次。前三次（代码等价、报告文字后续有修订）每次掉 1-2 个用例，全部落在
t1/t5/t6 已归类的「crash / interrupt / resume」flaky 家族；**第四次在最终提交上
一次全绿**：`514 passed, 25 skipped, 0 failed in 896.22s`（514+25 = 539 = 收集数）。

```
run 1 (892.22s): 1 failed, 512 passed, 25 skipped
  FAILED tests/test_adversarial_m3_challenger2.py::test_mid_run_crash_simulation_at_7_and_19_zero_loss_resume
run 2 (889.27s): 1 failed, 512 passed, 25 skipped
  FAILED tests/test_screen_runner.py::test_checkpoint_interrupt_simulation_and_zero_loss_resume
run 3 (892.72s): 2 failed, 512 passed, 25 skipped
  FAILED tests/test_adversarial_m3_challenger2.py::test_concurrent_crash_and_resumption
  FAILED tests/test_screen_runner.py::test_checkpoint_interrupt_simulation_and_zero_loss_resume
run 4 (896.22s, 最终提交 9397ecf): 0 failed, 514 passed, 25 skipped   ← 全绿
```

逐个复测（各 6 次，main 与本分支同样跑法）：

| 用例 | main（852d26c） | 本分支 |
|---|---|---|
| `test_mid_run_crash_simulation_at_7_and_19_zero_loss_resume` | 3 过 / 3 败 | 2 过 / 4 败 |
| `test_checkpoint_interrupt_simulation_and_zero_loss_resume` | 3 过 / 3 败 | 4 过 / 2 败 |
| `test_concurrent_crash_and_resumption` | 3 过 / 3 败 | 4 过 / 2 败 |

三个用例在 **main 上同样约 50% 失败** → 既有 flaky，与本次改动无关（与 t5 记录的
"base 2/6 败、head 3/6 败" 同量级）。三次全量运行的差异只体现在"这次轮到哪一个"。

- 新增测试文件 `server/tests/test_gnina_rescore_fix.py`：**508 行 / 16 个用例**；
  收集数：main `523` → run1/2 树 `538`（+15 新用例）→ 最终提交 `539`（+16）；每次
  都是 `passed + skipped + failed == collected`，无用例丢失。
- 相关子系统专项：`tests/test_screen_runner.py` + `tests/test_cli_screen.py` =
  **46 passed**（102.98 s）；`test_gnina_rescore_fix.py` + `test_gpu_adapters.py` +
  `test_scoring_rigor.py` = **51 passed / 2 skipped**（1.56 s）。

### 7.4 顺带修掉的一个"静默 skip"（T8-F2）

第一次按契约原样跑 verify-2（只设 `DOCKSTUDIO_REAL_GPU=1`，不导出 `DOCKSTUDIO_GNINA`）得到：

```
12 skipped, 1 warning in 0.06s
```

原因不是环境缺二进制，而是**测试自己的辅助函数找错了位置**：`_gnina_binary()` 只查
`<repo>/server/resources/bin/…`，而 `resources/bin` 其实在**仓库根**（`config.bundled_bin_dir()`
用 `parents[2]` 正是这个位置）。于是只要有 `DOCKSTUDIO_GNINA` 兜着就"看不见"这个问题，
一旦按契约写法跑就整套静默 skip —— 又一类"看着通过、其实没测"。

已修（`tests/test_gnina_real_gpu.py`：先查 `parents[2]` 再查 `parents[1]`），并把
`resources/bin/linux-x86_64/gnina` 复制进本工作树（gitignore 目录，不入库）。修后同一条
契约命令真的跑了：`10 passed, 2 xfailed in 32.21s`。

### 7.3 一个必须知道的环境陷阱（差点让本轮证据作废）

开发环境把 `dockstudio` 装成 **editable + meta-path finder**，指向 `/home/kui/dockstudio/server`。因此按契约写法 `cd /home/kui/wt/gnina-fix && pytest server/tests/…` 时：

```
$ cd /home/kui/wt/gnina-fix && python -c "import dockstudio; print(dockstudio.__file__)"
/home/kui/dockstudio/server/dockstudio/__init__.py        # ← 测的是 main！
```

第一次真机回归因此报了 2 failed（`DID NOT RAISE` + 老 `Parse error`），**看起来像修复没生效，其实是加载了 main 的代码**。这正是 t5 的 T5-F1。本轮在 `server/tests/conftest.py` 里把本工作树的 `server/` 插到 `sys.path[0]`（12 行），修复后：

```
DOCKSTUDIO_FILE: /home/kui/wt/gnina-fix/server/dockstudio/__init__.py
```

建议后续所有 worktree 的 pytest 证据都附这一行自证；`scripts/qa/run_gates.sh` 也应显式导出 `PYTHONPATH=<repo>/server`。

---

## 8. 遗留项与观察（不在本轮范围）

1. **GN-09 仍在**：gnina 的 `--minimize`/`--score_only` 输出里会出现**两个 `MODEL 1`**（本轮真机复现：`out.pdbqt` 第 1 行与第 6 行各一个 `MODEL 1`，只有一个 `ENDMDL`），解析后得到两个 `index=1` 的 pose。当前语义恰好正确（第一个是**优化后**的位姿，第二个是被回显的输入位姿），但 `Pose.index` 重复会让 UI/导出按 index 做 key 时冲突。属 t2 记录的 GN-09。
2. **GN-01 仍在**：配体文件不存在仍抛裸 `FileNotFoundError`（本轮未动，保持既有行为），真机测试里作为 `xfail` 保留。
3. **GN-05 仍在**：gnina 默认按 CNNscore 排序，`poses[0].affinity` 的语义与 vina 不同 —— 独立任务。
4. **GN-02/GN-06/GN-07/GN-10/GN-11/GN-12 未动**（stdout 五列表、报错截断、盒子前置校验、mock fixture、引擎名写死 vina、rmsd 恒空）。
5. **对 vina/smina/unidock 的行为变化（有意）**：这三个引擎现在也会做受体校验、也保留配体扩展名。前者只会把"空受体假成功"变成明确报错；后者只在调用方传入非 `.pdbqt` 配体时才与旧行为不同（旧行为必然解析失败）。

---

## 9. 复现步骤

```bash
# 0) 工作树
git -C /home/kui/dockstudio worktree add /home/kui/wt/gnina-fix -b fix/gnina-rescore main

# 1) 无 GPU 回归（快，1.8 s）
cd /home/kui/wt/gnina-fix/server && python -m pytest tests/test_gnina_rescore_fix.py -q

# 2) 门禁
cd /home/kui/wt/gnina-fix && ruff check server/dockstudio && ruff format --check server/dockstudio
cd /home/kui/wt/gnina-fix/server && python -m pytest tests/test_gpu_adapters.py tests/test_scoring_rigor.py -q

# 3) 真机（需 flock 独占 + gnina）
export DOCKSTUDIO_GNINA=/home/kui/dockstudio/resources/bin/linux-x86_64/gnina
flock /home/kui/.dockstudio-heavy.lock -c \
  "cd /home/kui/wt/gnina-fix && DOCKSTUDIO_REAL_GPU=1 python -m pytest server/tests/test_gnina_real_gpu.py -q -rxX"

# 4) Tier-3 端到端探针（4.1–4.4 的原始证据）
flock /home/kui/.dockstudio-heavy.lock -c \
  "DOCKSTUDIO_GNINA=$DOCKSTUDIO_GNINA python /home/kui/gnina_rescore_gpu_probe.py"
# 证据文件：/home/kui/gnina_probe_runs/rescore_fix_evidence.json
# 日志：/home/kui/gnina_probe_runs/rescore_fix_probe.log
```

> 说明：`/home/kui/gnina_rescore_gpu_probe.py` 是一次性探针（未入库），其断言已被 `server/tests/test_gnina_rescore_fix.py` 与 `test_gnina_real_gpu.py` 覆盖；证据 JSON/日志保留在机器上供复核。
