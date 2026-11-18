# 硬中止被静默吞掉 → 崩溃当正常完成（修复报告）

- 分支：`fix/abort-signal`（worktree `/home/kui/wt/abort`，基线 `main` = `a03beee`）
- 环境：kui@100.98.223.80，Python 3.11.15，pytest 9.2.0，ruff 0.16.10，Tesla P4
- 前置取证：`docs/qa/flaky-resume-2026-11-18.md`（t7，判定 (b) 产品侧缺陷，未改代码）
- 本报告：改了什么、为什么这样改、改造前后失败率、以及"崩溃/中止不再伪装成成功"的落点

---

## 0. 一句话结论

`ScreeningRunner` 曾用**同一个布尔** `_cancel_requested` 同时表示两件完全不同的事：
「调用方请收尾」（协作取消）和「进程要死了」（硬中止）。于是一次崩溃会走成这样：

```
worker: except BaseException → self._cancel_requested = True → re-raise
main  : ... → if self._cancel_requested: break   ← 先于下一次 fut.result()
        → 承载异常的 Future 从未被读取 → 异常被丢弃
        → run() 走到 _finalize_report() 正常返回一份截断报告
```

现在拆成两个信号：**硬中止一定抛出**（无论 Future 消费顺序），**协作取消仍是优雅返回**，
但返回的报告**不再看起来完整**（`funnel.csv` 多一行 `aborted`，且 T5 质量守恒不闭合）。

---

## 1. 改了什么（10 条 acceptance 逐条落点）

| # | 要求 | 落点 |
|---|---|---|
| 1 | 两个信号拆开，不再共用一个布尔 | 新增 `_hard_abort: BaseException \| None`、`_aborted_ids: set[str]`；`_cancel_requested` 只保留"协作取消"语义。`_record_hard_abort()` / `_raise_pending_abort()` / `_abort_gate()` 三个新方法 |
| 2 | 与 Future 消费顺序无关地抛出 | 机制见 §2；`_raise_pending_abort()` 在**每个**消费者循环之后、`_finalize_report()` 入口各调一次 |
| 3 | 覆盖单层 + 三层漏斗 | 单层 `_run_single_tier`（worker/consumer 都改）；三层 `_run_tier1_workers`、`_run_tier2_workers` 同构改造 |
| 4 | 两个 flaky 转确定性 | `test_abort_signal.py` 9 个新用例 + 两个既有用例实测 24/24 通过（失败率 100%→0%，见 §4） |
| 5 | 崩溃时报告不能看起来成功 | `funnel.csv` 新增 `aborted,<n>` 阶段行；`n_funnel_filtered` 不再把被丢弃的分子洗成"阈值淘汰"，因此 `is_conserved` 保持 False |
| 6 | 不回归 | `test_queue` / `test_job_lifecycle_e2e` / `test_htvs_funnel_stress` = 18 passed / 3 skipped（45.88s） |
| 7 | 证明加载本 worktree 代码 | `dockstudio.__file__ = /home/kui/wt/abort/server/dockstudio/__init__.py`（见 §5） |
| 8 | T9-F1（POSIX stub） | 两处 `#!/bin/sh` stub 改为 Python stub（`#!<sys.executable>`）；见 §6 |
| 9 | T9-F2（UniDock 批处理校验 + 过度声明） | `UniDockEngine.dock_batch` 接线 `validate_receptor_file()`（且**先于**二进制探测）；`base.py` docstring 改成准确表述 |
| 10 | T9-F7（gnina 报告数字） | `+1092/−42 → +1183/−49`、`54→53 files`、`a2ab77f→9397ecf`、补 `508 行 / 16 用例`、"四引擎全部继承"改成准确说法 |

改动文件（8 个）：

```
docs/fixes/gnina-rescore-2026-11-18.md    |  11 +-   （T9-F7 数字对齐 + T9-F2 措辞）
server/dockstudio/core/screen_runner.py   | 190 +++++-  （主修复）
server/dockstudio/engines/base.py         |  14 +-    （docstring 准确化）
server/dockstudio/engines/unidock.py      |  32 ++-    （批处理受体校验）
server/tests/test_abort_signal.py         | new       （9 个确定性回归用例）
server/tests/test_unidock_batch_guard.py  | new       （3 个批处理守卫用例）
server/tests/test_gnina_rescore_fix.py    |  83 ++-    （T9-F1 stub 跨平台）
server/tests/test_gpu_adapters.py         |  37 +-     （T9-F1 stub 跨平台）
```

---

## 2. 采用的机制，以及为什么它与完成顺序无关

三条互补的动作，缺一不可：

**(a) `except` 先存异常，再 re-raise（生产端）**

```python
try:
    res = self._dock_single_safe(lig, box, params, sub_dir)
except BaseException as exc:
    self._record_hard_abort(exc, lig.id)   # ← 第一件事：把异常本身记住
    self._cancel_requested = True
    raise
```

**(b) 消费者只对 `Exception` 降级，对 `BaseException` 记录后中止（消费端）**

```python
except (KeyboardInterrupt, SystemExit, GeneratorExit) as exc:
    self._record_hard_abort(exc, lig.id)
    self._cancel_requested = True
    self._mark_dock_failed(lig.id)
    executor.shutdown(wait=False, cancel_futures=True)   # 不再派发新任务
    break
except Exception as exc:            # 单个分子失败，运行继续
    ...
```

**(c) 循环之后 + finalize 入口无条件重放（顺序无关的关键）**

```python
        self._raise_pending_abort()      # ← 紧接 with 块之后

def _raise_pending_abort(self) -> None:
    exc = self._hard_abort
    if exc is None:
        return
    raise exc

# _finalize_report() 第一行再来一次（最后一道防线）
self._raise_pending_abort()
```

**为什么这样就与顺序无关**：`_hard_abort` 是 `self` 上的状态，**不依赖任何一个 Future 被读取**。
无论消费端先 `break`（原来的 bug 路径）、还是先读到那条 Future，
`_hard_abort` 都已被 worker 的 `except` 写过；(c) 的两处 `raise` 是**无条件**的，
所以"谁先谁后"不再影响结果。这就是 t7 报告 §6.1 的建议落地。

顺带修掉第二个问题（"崩了还白跑"）：worker 入口的 `_abort_gate()` 一旦发现
`_hard_abort` 已存在就立刻抛出，于是已提交但未开始的任务不再对接。
实测（12 个分子、第 3 个崩）：`dock_calls` 从**修前的 12** 降到**修后的 3**。

---

## 3. 协作取消 vs 硬中止：语义对照（可执行）

| | 协作取消 `runner.cancel()` | 硬中止（worker 抛 BaseException） |
|---|---|---|
| `run()` 返回 | **返回**部分报告（保持旧承诺） | **抛出**原异常 |
| 已完成的分子 | 保留在 checkpoint / 报告里 | 保留在 checkpoint / 报告里 |
| 未开始的分子 | 记入 `funnel.csv` 的 `aborted` 行 | 记入 `funnel.csv` 的 `aborted` 行 |
| `n_dock_failed` | 不变（不是引擎的错） | 崩溃的那个分子 +1 |
| `is_conserved` | **False**（不再伪装闭合） | 不适用（已抛异常） |
| CLI 表现 | `操作已取消 (User Cancelled)`，返回码 130 | `except KeyboardInterrupt` → 130，不打印"全部完成" |

实测（10 个分子、跑到第 4 个后取消）：

```
[single] docked_success=4  aborted=5  is_conserved=False  T5 expected=10 vs docked=4 -> OPEN
[tiered] docked_success=0  aborted=6  is_conserved=False  T5 expected=6  vs docked=0 -> OPEN
```

提交后（`99159bd`）的端到端复核，两条路径一起看：

```
[single] hard abort: RAISED KeyboardInterrupt  dock_calls= 4/10 funnel.csv=absent  checkpoint=3
[single] cancel    : RETURNED candidates=1     aborted_row=9  is_conserved=False
[tiered] hard abort: RAISED KeyboardInterrupt  dock_calls= 4/10 funnel.csv=absent  checkpoint=3
[tiered] cancel    : RETURNED candidates=0     aborted_row=9  is_conserved=False
```

两个额外的观察：

* **硬中止时 `funnel.csv` 根本不会被写**——`_finalize_report()` 在入口就重放了异常，
  所以磁盘上不会留下一份"看起来跑完了"的交付包。这比"写一份带 aborted 行的包"更干净。
* `dock_calls = 4/10`（10 个分子、第 4 个崩）——第 4 个之后一次都没再对接，
  印证了 `_abort_gate()` 确实止住了已排队的任务。

> 第二行是**修中间发现的一个真问题**：三层漏斗里被丢弃的分子原本会被算进
> `n_funnel_filtered`（"阈值淘汰"），于是 `expected == docked` 照样闭合、包看起来是完整的。
> 现已把 `aborted` 从 `funnel_filtered` 里扣除，T5 因此不闭合 —— 这才是正确行为。
>
> **但这套修法有漏洞，t12 审查抓到了（T12-F2），t14 已补**：上面那张表只对
> 「取消时还有分子**从未被 checkpoint**」成立。如果取消发生在**所有分子都已对接完**之后
> （tier 3 中途），`_abandoned_ligand_ids()` 返回 0，未评估的分子就被 `n_funnel_filtered`
> 吸收，`is_conserved` 仍然是 True —— 详见 §11.2（t14 追加）。

---

## 4. 失败率：改造前 → 改造后（各 ≥20 次重复）

跑法：每次一个独立进程，`PYTHONPATH=<本 worktree>/server`，**不加 sleep、不改断言、不 deselect**。

| 用例 | 改造前（24 次） | 改造后（24 次） |
|---|---|---|
| `test_adversarial_m3_challenger2.py::test_mid_run_crash_simulation_at_7_and_19_zero_loss_resume` | **18 败 / 6 过 = 75.0% 败** | **0 败 / 24 过 = 0.0%** |
| `test_screen_runner.py::test_checkpoint_interrupt_simulation_and_zero_loss_resume` | **11 败 / 13 过 = 45.8% 败** | **0 败 / 24 过 = 0.0%** |
| 两个文件一起（verify 命令 2 的口径） | — | **24/24 全绿** |

改造前的失败首行一律是 `Failed: DID NOT RAISE KeyboardInterrupt`
（`test_adversarial_m3_challenger2.py:207` / `:227`、`test_screen_runner.py:674`），
即本报告要修的"濒死进程返回了正常报告"。

新增的 9 个用例（`server/tests/test_abort_signal.py`）**不靠时序**：
mock 引擎在 `dock()` 里用 `threading.Event` 打闸，由测试决定放行时刻，
所以"消费端先看到 flag"和"Future 从未被读"两种顺序都被**强制**跑过一遍：

```
test_hard_abort_propagates_single_tier
test_hard_abort_propagates_when_consumer_already_saw_cancel
test_run_resets_abort_state_between_runs
test_hard_abort_propagates_tier1_tiered_funnel
test_hard_abort_does_not_drain_worker_queue
test_cooperative_cancel_returns_partial_report_but_marks_truncation
test_cooperative_cancel_is_not_reported_as_dock_failure
test_abort_then_resume_zero_loss_three_phases
test_tiered_abort_then_resume_reaches_full_completion
```

---

## 5. 测试加载的是本 worktree 的代码（T5-F1 那个坑）

```
$ PYTHONPATH=/home/kui/wt/abort/server python -c "import dockstudio; print(dockstudio.__file__)"
/home/kui/wt/abort/server/dockstudio/__init__.py
```

本分支的 `server/tests/conftest.py`（t8 加的）会把 `<worktree>/server` 插到 `sys.path[0]`，
所以即使忘了 `PYTHONPATH` 也会加载本树；上面这行是在**显式** `PYTHONPATH` 下取的证。
本次所有 pytest 命令都带 `PYTHONPATH=/home/kui/wt/abort/server`。

---

## 6. T9-F1：POSIX shell stub → 跨平台 Python stub

**具体哪几处依赖 POSIX shell**（改前）：

| 位置 | 依赖 |
|---|---|
| `server/tests/test_gnina_rescore_fix.py:112` `_fake_gnina()` | `#!/bin/sh` + `for a in "$@"` + `grep` + `printf` + `case … esac` |
| `server/tests/test_gpu_adapters.py:136` `test_gnina_mock_binary_e2e` | `#!/bin/sh` + `for a in "$@"` + `printf` |
| 两处 `script.chmod(0o755)` | 依赖可执行位 + shebang（Windows 无 `sh`，`chmod` 也近乎无效） |

`ci.yml` 的 `server-test-full` 在 `windows-latest` 上跑同一批测试，所以这两处会让全矩阵新增失败。

**修法**：把两处 stub 改成 **Python 脚本**，shebang 用 `#!<sys.executable>`（不是 `env python`：
某些环境只有 `python3` 在 PATH 上，Windows 更没有 `env`）。行为完全保留：
`--version` 短接、`argv.txt` 落盘、PDBQT 遇到 `V2000` / `ENDMDL` 报同样的两行 Parse error、
输出 `MODEL 1` + CNN REMARK。`chmod` 保留（POSIX 上仍需要可执行位）。
顺带删掉 `test_gpu_adapters` 里那段 `if platform.system() == "Windows" or True: pass` 死代码。

改后：`test_gnina_rescore_fix.py` + `test_gpu_adapters.py` = **25 passed / 2 skipped**。

> ⚠️ **本节原标题写的「windows-latest 不再新增失败」是错的，t12 审查已判为 T12-F1，此处更正。**
>
> 事实：Python 脚本的 `#!<python>` shebang 只在 POSIX 上由内核解释。Windows 的
> `CreateProcess` **不经过 shell、也不读 shebang**，所以把一个 `.py` 文件当可执行程序
> 直接 `subprocess` 起来会失败 —— 换成 Python stub 只解决了「不需要 `sh`」，
> 没有解决「Windows 上这个文件根本启动不了」。
>
> 判定依据（实测，不是推理）：在 Linux 上把 stub 的 shebang 删掉，内核同样返回
> `ENOEXEC`（`OSError: [Errno 8] Exec format error`），这与 Windows 上的失败模式一致。
> 用 pytest 插件在测试期做这个替换后跑 `test_gnina_rescore_fix.py`：
> **5 failed / 11 passed**，失败的正好是下面这 5 个 ——
> `test_in_place_rescore_accepts_an_exported_sdf_pose`、
> `test_multimodel_out_pdbqt_is_reduced_to_one_model`、
> `test_empty_receptor_is_rejected_before_launching`、
> `test_receptor_with_atoms_is_accepted`、
> `test_cnn_missing_result_carries_a_degradation_marker`
> （另有 `test_gpu_adapters.py::test_gnina_mock_binary_e2e` 同样依赖可执行 stub）。
>
> 修法：给这 6 个用例加 `@pytest.mark.skipif(sys.platform == "win32", ...)`，理由写在
> skip reason 里。**不**提供 Windows launcher，因为那要么引入 `.bat`（又一次把测试绑到
> 平台 shell 上），要么得给引擎适配器加一个"用 `sys.executable` 起脚本"的后门 ——
> 而那个后门只在测试里有意义，会把生产路径搞脏。这些用例覆盖的是 Linux/macOS 上的
> gnina 交互，Windows 侧本来就没有真机跑道（CI 的 windows-latest 只跑零 GPU 单测）。
>
> 验证：模拟 `sys.platform == "win32"` 时该文件 **19 passed / 8 skipped**（其中 6 个是
> 本 finding 的 skip，2 个是既有的 off-linux skip），Linux 上仍是 **25 passed / 2 skipped**。

---

## 7. T9-F2：UniDock 批处理入口的受体校验 + 去掉过度声明

**改前的事实**：`validate_receptor_file()` 只在 `VinaEngine.dock()` 里被调用。
`UniDockEngine.dock_batch()` 是**独立入口**，不经过 `dock()`；而 `base.py` 的 docstring 写
"Both the generic adapters and the HTVS runners call this"，报告里写"四个引擎全部继承到" ——
两句都不准确。

**改法**：

```python
# unidock.py dock_batch()，注意顺序
# 1. Argument-order polymorphism support（纯解析，无子进程）
# 2. validate_receptor_file(receptor_path)   ← 新增，先于二进制探测
# 3. Check engine availability（_binary()/_version()）
```

> 顺序很重要：如果先探测二进制，一台**没装 Uni-Dock** 的机器会先抛
> `EngineUnavailableError`，把"受体是空的"这个真正的输入错误盖掉，而且这个守卫就只在
> 装了 GPU 引擎的机器上才有效 —— 对 CI 恰好是反的。

`base.py` docstring 改成列出**真实**调用点（`VinaEngine.dock`、`UniDockEngine.dock_batch`），
并明说 `BaseEngine.dock`/`dock_batch` 是协议/兜底、**不会**自动调用（新引擎必须显式接线）。

新增 `server/tests/test_unidock_batch_guard.py`（3 个用例）：空受体 → `DockingError`（含"原子"）、
不存在的受体 → `DockingError`、有效受体 → 通过守卫（不再报"原子"）。

---

## 8. T9-F7：gnina 报告的与事实不符的数字

| 位置 | 改前 | 改后 | 依据 |
|---|---|---|---|
| 头部 diffstat | `9 文件 +1092 / −42` | `9 文件 +1183 / −49` | `git show --stat 9397ecf` |
| §7.1 ruff | `54 files already formatted` | `53 files already formatted` | 本树 ruff format --check |
| §7.2 全绿那次的 commit | `a2ab77f` | `9397ecf` | 分支 tip（已合入 main） |
| §7.2 测试文件规模 | 未写 | `508 行 / 16 个用例` | `git show 9397ecf:server/tests/test_gnina_rescore_fix.py` |
| §1 GN-04 行 | "四个引擎全部继承到" | 准确表述 + T9-F2 备注 | 见 §7 |

---

## 8.1 全量 pytest 结果（两种条件，都不回归）

| 条件 | 结果 |
|---|---|
| 本 worktree，**未注入引擎**（干净工作树无 `resources/bin`） | `1 failed, 513 passed, 37 skipped in 307.09s` —— 唯一失败 `test_cli_agent.py::test_dock_demo_and_post_dock_actions_end_to_end`，首行 `EngineUnavailableError: vina binary not found; set DOCKSTUDIO_VINA or bundle vina under resources/bin`（**环境产物**，与本次改动无关；t1 已记录同一现象） |
| 注入 `DOCKSTUDIO_VINA`（= **CI 真实条件**，`fetch_engines.sh` 会放好二进制） | **`526 passed, 25 skipped, 0 failed in 921.20s`（exit 0）（99159bd 上复跑：526 passed / 25 skipped / 0 failed / 899.34s / exit 0；t14 修完 T12-F1/F2 后为 531 passed / 25 skipped / 0 failed / 962.76s，见 §11）** |

25 个 skip 的构成全部是既有 opt-in / 平台守卫，无一新增：

```
12  test_gnina_real_gpu.py     real-GPU test is opt-in (DOCKSTUDIO_REAL_GPU=1)
10  tests/golden/test_golden.py golden set is opt-in (-m golden)
 2  test_gpu_adapters.py       windows/off-linux check
 1  test_smina.py              smina binary not available
```

> 收集数：main `539` → 本分支 `551`（+9 `test_abort_signal.py` + 3 `test_unidock_batch_guard.py`）。
> 两次运行都满足 `passed + skipped + failed == collected`，无用例丢失。
> **注意**：本次没有用 `DOCKSTUDIO_GATES_DESELECT`，也没有环境变量包装 flaky 用例 ——
> 修好之后它们在两种条件下都自然通过。

---

## 9. 复现命令（三条 verify 命令 + 关键取证）

```bash
PY=/home/kui/miniconda3/envs/ai-gpu/bin/python
WT=/home/kui/wt/abort
LOCK=/home/kui/.dockstudio-heavy.lock

# ① 不回归三件套
flock $LOCK -c "cd $WT && PYTHONPATH=$WT/server $PY -m pytest \
  server/tests/test_queue.py server/tests/test_job_lifecycle_e2e.py \
  server/tests/test_htvs_funnel_stress.py -q"
# 18 passed, 3 skipped in 45.88s  (exit 0)

# ② 24 次重复（两个曾经 flaky 的文件）
flock $LOCK -c 'cd /home/kui/wt/abort && for i in $(seq 1 24); do \
  PYTHONPATH=/home/kui/wt/abort/server '"$PY"' -m pytest \
  server/tests/test_adversarial_m3_challenger2.py server/tests/test_screen_runner.py -q \
  || echo FAIL_RUN_$i; done'
# 24/24 全绿

# ③ 静态门禁
cd $WT && /home/kui/miniconda3/envs/ai-gpu/bin/ruff check server/dockstudio     # All checks passed!
cd $WT && /home/kui/miniconda3/envs/ai-gpu/bin/ruff format --check server/dockstudio  # 53 files already formatted

# 额外：新回归套件 / 取消语义 / 排空队列
flock $LOCK -c "cd $WT && PYTHONPATH=$WT/server $PY -m pytest server/tests/test_abort_signal.py -q"        # 9 passed
flock $LOCK -c "cd $WT && PYTHONPATH=$WT/server $PY -m pytest server/tests/test_unidock_batch_guard.py -q" # 3 passed
```

---

## 10. 边界与未做的事

- **只改了产品侧的三个文件**（`core/screen_runner.py`、`engines/base.py`、`engines/unidock.py`），
  其余全是测试与文档；`apps/`、`resources/`、`scripts/`、`docs/qa/` 未动。
- **没有跑真机 GPU 用例**（`DOCKSTUDIO_REAL_GPU=1`，12 条按设计 skip）：本次改动不触 GPU 路径，
  但也没有验证过 GPU 批处理分支上的中止行为。
- **`_run_single_tier_batch` / `_run_tier1_batch` / `_run_tier2_batch` 未做硬中止改造**：
  它们在主线程串行调用 `engine.dock_batch()`，一旦抛 `BaseException` 本来就会直接穿出 `run()`，
  不存在"异常落进 Future 没人读"的窗口；但这也意味着"批处理中途崩溃"不会把已完成的
  chunk 记账（另一个方向的问题，本轮未处理）。
- **协作取消的对外可见性**只做到 `funnel.csv` 的 `aborted` 行 + `is_conserved`
  （t14 把这条补全到"每个 tier 都认账"，见 §11.2）；
  没有新增 `ScreeningFunnelReport.aborted` 字段（那会改 schema，属于需要单独评审的 API 变更）。
  CLI 侧沿用既有的 `KeyboardInterrupt → 130` 分支，未改 UI 文案。
- **Windows 侧没有真机跑道**：T12-F1 的做法是 skip（理由见 §11.1），不是让那 6 个用例
  在 Windows 上跑起来。如果以后要给 Windows 加 gnina 交互测试，需要一个真正的
  launcher 方案（而不是 `.py` + shebang）。
- **`_run_*_batch` 的截断记账仍是空白**：GPU 批处理入口一次吃一个 chunk，chunk 中途崩
  会丢掉该 chunk 已完成的分子（本轮只覆盖了线程池路径）。
- **改造前的失败率随机器负载浮动**：本机实测 75.0% / 45.8%（t7 记录 62.5% / 45.8%）；
  改造后是 0/24 而不是"1-x"，判据应该看 `slow-main` 那种强制时序的场景（新用例已覆盖）。

---

## 11. t14 追加：补齐 t12 审查的两条 medium

t12 复核判 pass，但留了 4 条 findings；t14 把其中两条 medium 与两条 low 全部处理掉。
本节是增量记录，前面的 §1–§10 保持原样（除 §6 的更正与 §8.1 的编号）。

### 11.1 T12-F1：Windows 上那 5 个用例其实会红 —— 更正 §6 的结论

§6 原写「windows-latest 不再新增失败」，**这句是错的**，已在 §6 内就地更正并给出判定依据。
一句话复述：Python 脚本的 shebang 只在 POSIX 由内核解释，Windows 的 `CreateProcess`
不经过 shell、不读 shebang，所以 `.py` stub 在 Windows 上根本无法作为程序启动。

- 判定依据（实测）：Linux 上删掉 shebang → 内核返回 `ENOEXEC`（`[Errno 8] Exec format error`），
  与 Windows 失败模式一致；用 pytest 插件在测试期做这个替换 → `test_gnina_rescore_fix.py`
  **5 failed / 11 passed**，失败集合与 T12-F1 点名的完全一致。
- 修法：6 个用例（含 `test_gpu_adapters.py::test_gnina_mock_binary_e2e`）加
  `@pytest.mark.skipif(sys.platform == "win32", ...)`，理由写进 skip reason。
- 为什么不提供 Windows launcher：`.bat` 等于再绑一次平台 shell；给适配器加"用
  `sys.executable` 起脚本"的后门只为测试存在，会污染生产路径。这些用例验证的是
  gnina 进程交互，Windows 侧本来没有真机跑道。
- 验证：模拟 `sys.platform == "win32"` → 该两文件 **19 passed / 8 skipped**（6 条本 finding +
  2 条既有 off-linux）；Linux 上仍是 **25 passed / 2 skipped**。

### 11.2 T12-F2：协作取消也会"看起来守恒"—— 已修，且有数字

**t12 报的场景我先复现了，并确认它比报告的还宽**（n=12，`tier1_keep_ratio=tier2_keep_ratio=1.0`）：

| 场景 | 修前 candidates | 修前 funnel_filtered | 修前 is_conserved | 修前 aborted 行 |
|---|---|---|---|---|
| tier 3 中途取消 | 1 | **11** | **True** ❌ | 无 |
| tier 2 中途取消（keep1=0.5） | 0 | **12** | **True** ❌ | 无 |
| tier 3 中途取消（keep1=0.5） | 1 | **11** | **True** ❌ | 无 |

根因不是"取消没记账"，而是**记账的判据用错了层级**：`_abandoned_ligand_ids()` 问的是
"这个分子有没有**任何** checkpoint 记录"。取消发生在 tier 2/3 时，那些分子早在 tier 1
就 checkpoint 过了，于是被判为"处理过"，未评估的部分被 `n_funnel_filtered` 吸收
（它本来是"阈值淘汰"的桶），守恒式照旧闭合。

**修法（选"强制不守恒 + 写 aborted 行"，不是收窄措辞）**：
按 tier 记录截断，而不是按"有没有 checkpoint"：

```python
self._tier_truncated: dict[int, set[str]] = {1: set(), 2: set(), 3: set()}

def _record_tier_truncation(self, tier, pending_ids, completed_ids):
    if not self._cancel_requested:
        return
    missing = {str(i) for i in pending_ids if str(i) not in completed_ids}
    self._tier_truncated.setdefault(tier, set()).update(missing)
```

三个 tier 的 worker 循环结束后各调一次；`run()` 的两个早返回分支（tier1 后 / tier2 后）
也会把"上一层的幸存者"记为下一层的截断 —— 否则取消正好落在 tier 边界时，
那批分子谁都不认领（实测会少报 1 个：`aborted=5` 而实际 6 个未评估）。

记账侧：`n_tier1_filtered`/`n_tier2_filtered` 各自减去本层截断数，
`n_funnel_filtered = tier1_filtered + tier2_filtered`（不再用 `pool - success - failed` 的减法），
截断集合与 `_abandoned_ligand_ids()` 合并后写 `aborted` 行。

**为什么不用"残差兜底"**：第一版我用 `residual = pool - success - failed - filtered` 把
无法解释的分子一律算成截断。这在完整运行上会**误伤** —— tier 3 的 QC 淘汰
（如 `HIGH_STRAIN_UNCONVERGED_FALLBACK`）本来就不计入 `funnel_filtered`，
残差会把它们标成"被中止"，让一次正常跑完的筛选手动变成"看起来被取消"。
所以最终版**只认能证明的截断**（进了这一层、没出来），不做兜底猜测。

修后（同一组场景）：

| 场景 | candidates | funnel_filtered | is_conserved | aborted 行 |
|---|---|---|---|---|
| tier 3 中途取消 | 1 | 0 | **False** ✅ | **11** |
| tier 2 中途取消（keep1=0.5） | 0 | 6（= 真实 12→6 比例淘汰） | **False** ✅ | **6** |
| tier 3 中途取消（keep1=0.5） | 1 | 6 | **False** ✅ | **5** |
| **完整跑完**（不取消） | 12 | 0 | **True** ✅ | 无 |

最后一行是关键的反向对照：不取消时守恒式仍然闭合、不写 `aborted` 行 ——
修法没有把"正常完成"误判成"被截断"。

新增 4 个回归用例（`test_abort_signal.py`）：
`test_cancel_during_tier3_does_not_claim_conservation`、
`test_cancel_during_tier2_does_not_claim_conservation`、
`test_cancel_accounting_balances`（`filtered + failed + success + aborted == raw`）、
`test_complete_run_still_conserves`。

### 11.3 T12-F3：硬中止不再留下上一次的交付包（两行改动）

复用 `--out` 时，硬中止会让磁盘上留着**上一次完整运行**的 `funnel.csv` /
`screening_results.csv` —— 一个中止的目录里躺着一份"看起来完成了"的交付包。
修法：`_record_hard_abort()` 里顺手删掉这两个文件（`unlink(missing_ok=True)`，
失败只记 warning，不让清理动作压过中止本身）。新增用例
`test_hard_abort_removes_stale_delivery_files`。

### 11.4 T12-F4：报告自身的三处不实

| 项 | 改前 | 改后 |
|---|---|---|
| 行数自述 | 326 行（实际 325） | 本版 **480 行**（`wc -l` 实测） |
| 章节编号 | `§8` 后接 `§8.5`（重复层级） | `§8.1` |
| commit 引用 | `9accc37`（amend 前的旧 hash） | `99159bd` |

### 11.5 t14 的验收数字

| 项 | 结果 |
|---|---|
| `test_abort_signal.py` + `test_screen_runner.py` + `test_adversarial_m3_challenger2.py` | **58 passed**（10.60s） |
| 全量 pytest（注入 `DOCKSTUDIO_VINA` = CI 真实条件） | **531 passed / 25 skipped / 0 failed / 962.76s / exit 0** |
| 两条原 flaky 各 24 次隔离运行 | **24/24 + 24/24，失败率 0.0%** |
| `ruff check server/dockstudio` / `ruff format --check` | `All checks passed!` / `53 files already formatted` |
| 新增用例 | `test_abort_signal.py` 9 → **14** 个 |

> 收集数从 551 增至 556（+5 个新用例：4 个 T12-F2、1 个 T12-F3），
> `passed + skipped + failed == collected` 成立，无用例丢失。
> 仍未使用 `DOCKSTUDIO_GATES_DESELECT`，也没有环境变量包装。

---
