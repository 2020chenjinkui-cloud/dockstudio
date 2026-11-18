# 崩溃续跑用例不稳定：根因定位报告（2026-11-18）

- 被测用例：`tests/test_adversarial_m3_challenger2.py::test_mid_run_crash_simulation_at_7_and_19_zero_loss_resume`
- 基线 commit：`c0f33cc7dbea2928d9840a18f4ab11c66d8bc317`（main）
- 环境：kui@100.98.223.80，Python 3.11.15，pytest 9.2.0，Linux 7.0.0-38-generic
- 取证脚本：`scripts/qa/repro_flaky_resume.py`（本次新增，只读产品代码）
- **本轮未修改 `server/dockstudio/**` 任何一行**；本报告只定位与取证

---

## 0. 结论（先给判定）

**判定：(b) 产品侧缺陷。** 但必须把"坏在哪"说准：

| 承诺面 | 结论 | 证据 |
|---|---|---|
| **checkpoint 数据完整性 / 续跑零丢失** | **成立**（未发现丢分子、重复、半行） | 53 次逐次打印不变量的迭代里，`dock_calls` 恒为 7、`checkpoint_lines` 恒为 6 |
| **中止信号传播（异常）** | **不成立 —— 会被静默吞掉** | 62.5%（15/24）的隔离运行里 `run()` 正常返回，`KeyboardInterrupt` 从未被消费 |

一句话机制：

> 单层对接路径里，worker 在自己的 `except BaseException` 里把 `self._cancel_requested` 置 True **之后**才 re-raise；
> 而消费端循环末尾的 `if self._cancel_requested: break`（`screen_runner.py:1524`）**先于**下一次 `fut.result()` 生效。
> 一旦 worker 抢在主线程读到那个"成功结果"之前/之际置位，主线程就 `break` 走人，**承载崩溃异常的那个 Future 从头到尾没有被取过**，
> 异常随对象一起被丢弃，`run()` 平稳走到 `_finalize_report()` 返回一份截断报告。

这不是"测试发信号的时机不对"：用例里没有任何 sleep/时序等待，mock 引擎的崩溃注入是 `raise KeyboardInterrupt(...)` 的**确定性**动作（`test_adversarial_m3_challenger2.py:95`）。竞态变量完全在产品代码里——**worker 写 flag 与主线程读 flag 的先后**。

---

## 1. 失败率（实测数字）

### 1.1 pytest 隔离跑该用例（每次一个新进程）

```bash
cd /home/kui/wt/qa-baseline/server
for i in $(seq 1 24); do
  /home/kui/miniconda3/envs/ai-gpu/bin/python -m pytest \
    tests/test_adversarial_m3_challenger2.py::test_mid_run_crash_simulation_at_7_and_19_zero_loss_resume \
    -q --tb=line
done
```

| 项 | 数字 |
|---|---|
| 运行次数 N | **24** |
| 失败 | **15** |
| 通过 | **9** |
| **失败率** | **62.5%** |
| 失败点 A：`tests/...py:207` → `Failed: DID NOT RAISE KeyboardInterrupt`（Run 1 在 LIG_007 崩溃被吞） | **9 次** |
| 失败点 B：`tests/...py:227` → 同上（Run 2 在 LIG_019 崩溃被吞） | **6 次** |
| 数据断言类失败（`len(lines)==6/18/30`、重复、candidates 数量） | **0 次** |

> 失败点 B 是本轮新发现：不只 Run 1 会吞，**续跑后的 Run 2 同样会吞**，说明缺陷在代码路径上而不是某个一次性初始化。
> 「0 次数据断言失败」也是硬证据：**只要异常能抛出来，所有 checkpoint/零丢失断言都通过**。

历史观测（前序任务 t1）：全量套件里 7 次混合运行 = 1 胜 6 败（全量 Run A 通过、Run B 失败、单独跑 5 次全败）。

#### 1.1.1 失败堆栈摘录（原样）

失败点 A（Run 1 在 LIG_007 的崩溃被吞）：

```
>       with pytest.raises(KeyboardInterrupt, match="Simulated abrupt process crash at ligand LIG_007"):
E       Failed: DID NOT RAISE KeyboardInterrupt
tests/test_adversarial_m3_challenger2.py:207: Failed
=========================== short test summary info ============================
FAILED tests/test_adversarial_m3_challenger2.py::test_mid_run_crash_simulation_at_7_and_19_zero_loss_resume
1 failed, 1 warning in 1.61s
```

失败点 B（续跑后的 Run 2 在 LIG_019 的崩溃同样被吞）：

```
>       with pytest.raises(KeyboardInterrupt, match="Simulated abrupt process crash at ligand LIG_019"):
E       Failed: DID NOT RAISE KeyboardInterrupt
tests/test_adversarial_m3_challenger2.py:227: Failed
```

同一次运行的被测代码日志（`-o log_cli=true --log-cli-level=INFO`）——注意最后一行：**6 个候选、Conserved: False，但没有任何异常**：

```
INFO  dockstudio.core.screen_runner: Starting screening run in .../crash_resume_run
      (engine=mock_crash_engine, resume=False, tiered=False)
INFO  dockstudio.core.screen_runner: Single-tier: 0 completed from checkpoint, 30 pending
WARNING posebusters.modules.distance_geometry: Molecule has only 1 atoms.   # x6（LIG_001..006 各一次）
INFO  dockstudio.core.screen_runner: Screening report finalized: 6 final candidates.
      Conserved: False. Output: .../crash_resume_run
```

正常抛出时（同一用例的 PASS 运行）对照堆栈：

```
KeyboardInterrupt: Simulated abrupt process crash at ligand LIG_007
  File ".../dockstudio/core/screen_runner.py", line 1397, in run
    return self._run_single_tier(
  File ".../dockstudio/core/screen_runner.py", line 1516, in _run_single_tier
    res_dict = fut.result()
  File ".../dockstudio/core/screen_runner.py", line 1505, in _do_single_dock
    res = self._dock_single_safe(lig, box, params, sub_dir)
  File ".../dockstudio/core/screen_runner.py", line 2470, in _dock_single_safe
    return self.engine.dock(
  File ".../tests/test_adversarial_m3_challenger2.py", line 95, in dock
    raise KeyboardInterrupt(f"Simulated abrupt process crash at ligand {lig_id}")
```

### 1.2 取证脚本逐档对比（`scripts/qa/repro_flaky_resume.py`）

| 模式 | 做法 | 结果 |
|---|---|---|
| `bare` | 完全不动时序 | 20 次中 **7 次 NO_EXCEPTION（35%）**、13 次正常抛出 |
| `slow-main` | 主线程在 `_emit_progress` 里 sleep 0.3s（`_emit_progress` 就是 `1524` 那个检查点前最后一个动作）→ **放大「flag 已置位」窗口** | **6/6 全部 NO_EXCEPTION**（另有早前 probe 3/3，合计 9/9） |
| `slow-worker` | worker 每次 `_dock_single_safe` 前 sleep 0.3s → **把 worker 拖到消费端后面** | **6/6 全部正常抛出**（另有早前 probe 3/3，合计 9/9） |
| `boundary` | 边界打点（见 §3） | 见 §3 |

`slow-main` 100% 吞没 / `slow-worker` 0% 吞没 —— 同一个缺陷，两个方向一拧就翻转，**这是"消费端先于生产端观测 flag"的直接证明**。

---

## 2. 用例期望的时序 vs 实际时序

```
用例期望（单层路径，workers = max(1, cpu//cpu_per_job) = max(1, 2//2) = 1）：
  worker: dock(LIG_001..006) → 各自 checkpoint.record_single()
          dock(LIG_007) → raise KeyboardInterrupt("Simulated abrupt process crash at ligand LIG_007")
  main  : fut.result() 依次取 001..006 → 每条都 _emit_progress
          fut.result() 取 007 → 抛出 KeyboardInterrupt → pytest.raises 接住 → PASS

实际（约 62.5% 的运行）：
  worker: ... dock(LIG_006) 完成 → 立刻开始 dock(LIG_007) → 崩溃 → self._cancel_requested = True（1507 行）
  main  : 还在处理 LIG_006 的返回值 / 刚进 _emit_progress
          → 走到 1524 行 `if self._cancel_requested:` → True → break
          → 循环结束，007 的 Future 从未被 result() 取过 → 异常被丢弃
          → 继续 _finalize_report() → run() **正常返回**（candidates=6，is_conserved=False）
```

要紧的是：`_do_single_dock` 的 `except BaseException`（1506-1508）把"硬崩溃"和"协作取消"混成了同一个布尔标志 `_cancel_requested`，而消费端把这个标志当作"别再看结果了，收摊"的信号。

---

## 3. 证据链（5 条，逐条可复现）

### E1 — 不改任何东西就能复现"静默返回"

`--mode bare`：20 次里 7 次 `run()` 正常返回（`NO_EXCEPTION`），13 次抛 `KeyboardInterrupt`。
每次返回的报告都是：`n_ready=30, n_docked_success=6, n_dock_failed=0, candidates=6, is_conserved=False`。

→ 崩溃点之后的 23 个分子根本没被处理，但 `run()` 没有任何异常。**调用方无法从返回值/异常区分"崩溃中止"与"跑完了"。**

### E2 — 放大消费端窗口 → 吞没率 100%

`--mode slow-main --delay 0.3`：9/9 全部 `NO_EXCEPTION`（其中一次 `candidates` 只剩 **2**，而 checkpoint 里明明是 6 条）。
→ 说明窗口就在"worker 置位"和"主线程检查"之间；把消费者拖慢，吞没变成必然而不是偶然。

### E3 — 反向推迟生产端 → 吞没率 0%

`--mode slow-worker --delay 0.3`：9/9 全部抛出 `KeyboardInterrupt`。
→ 同一份产品代码、同一个用例，只改线程时序就翻转结果 —— 竞态定位成立。

### E4 — 边界打点：异常"从未被消费"（最硬的证据）

对 `Future.result` / `ScreeningRunner._mark_dock_failed` / `_run_single_tier` / `run` 四个边界打点：

吞没场景的事件序列（完整）：

```
_run_single_tier RETURNED in MainThread
run RETURNED in MainThread
```

抛出场景的事件序列（完整）：

```
Future.result RAISED KeyboardInterrupt in MainThread
_mark_dock_failed(LIG_007) in MainThread
_run_single_tier RAISED KeyboardInterrupt in MainThread
run RAISED KeyboardInterrupt in MainThread
```

即：吞没时 **`Future.result` 一次都没抛、`_mark_dock_failed(LIG_007)` 一次都没被调用** ——
`1518` 行那个 `except BaseException as exc:` 分支压根没进过；循环是通过 `1524` 行的 `break` 退出的。

控制流落点（吞没场景下 `_finalize_report` 被调用时的栈）：

```
File ".../dockstudio/core/screen_runner.py", line 1397, in run
    return self._run_single_tier(
File ".../dockstudio/core/screen_runner.py", line 1614, in _run_single_tier
    return self._finalize_report(
```

→ 确认是从 `with ThreadPoolExecutor(...)` 块之后**正常落下**到 1614 行，而不是异常穿出。

### E5 — 对照：三层漏斗路径不会吞（缺陷只在单层路径）

同样的 mock 崩溃引擎，但 `enable_tiered_funnel=True`：**5/5 全部正常抛出**。
代码层面也自洽：`_run_tier1_workers` 的 `_do_dock`（1763-1791）不捕获 `BaseException`、不置位 flag，
消费端 `1795` 的 `if self._cancel_requested: break` 因此不会触发，`fut.result()` 的 `KeyboardInterrupt` 直接穿出。

> 顺带发现（同一份代码，另一个小问题）：三层路径下崩溃后**队列不排空**——worker 会把剩下所有分子都 dock 完再抛
> （探针里 `dock_calls=10`，10 个分子全跑完才抛）。这是"崩了还白跑"，与本次"吞异常"是两件事，仅记录。

---

## 4. 为什么判 (b) 而不是 (a)（测试侧竞态）

| 判据 | 事实 |
|---|---|
| 用例里有没有"保守等待/sleep"可调？ | **没有**。用例零 sleep，mock 引擎的崩溃是确定性 `raise`（`:95`） |
| 断言是不是产品对外契约？ | 是。`test` 断言的正是"崩溃 → `run()` 必须抛 `KeyboardInterrupt`"，同一份意图也写在 `screen_runner.py:14` 的模块 docstring（"Zero-loss continuation on restart with resume=True"）里 |
| 失败是否伴随数据错误？ | **从不**。15 次失败全是 `DID NOT RAISE`，`len(lines)` / 重复 / candidates 断言 0 次失败 |
| 竞态变量在哪？ | 产品代码：worker 写的 `_cancel_requested`（1507）与主线程读的 `_cancel_requested`（1524），两者之间没有 happens-before 约束 |
| 反向验证 | `slow-main` 100% 吞、`slow-worker` 0% 吞（E2/E3）—— **改的是"哪一侧"决定结果，产品行为不是测试构造出来的** |

因此：**调 sleep / 放宽断言 / 删用例都不是修复**（也符合本任务"不得为了让测试变绿而调 sleep"的约束）。
改测试唯一能做到的是"让它永远不测这条路径"，等于把缺陷藏起来。

---

## 5. 对"断电续跑零丢失"这一对外承诺的影响

必须分开说，因为两条腿一好一坏：

**好的那条腿（数据）：成立。**
- 53 次逐次打印不变量的迭代（含 ≥18 次吞没场景）里：`dock_calls` 恒为 **7**（LIG_001..007 各 1 次，崩溃后无重试、无重复）、`checkpoint.jsonl` 恒为 **6** 行（LIG_001..006）。
- 一旦异常正常抛出，用例的 6 / 18 / 30 三档断言全部通过 ⇒ **续跑不会重复对接、不会丢已完成分子**。
- 也就是说："进程死掉 → checkpoint 里已完成的就是已完成的 → 下次 resume 从这里继续" 这条承诺**代码是守的**。

**坏的那条腿（信号）：不成立。**
- 崩溃中止可能被静默吞掉，`run()` 返回一份 `candidates=6`、`is_conserved=False` 的报告。
- 靠异常判中止的调用方（CLI 的 `except KeyboardInterrupt` 分支、agent、队列封装、CI 脚本）会把它当作**正常完成**。
- 更糟的次生不一致：吞没时 `n_dock_failed=0`（LIG_007 连"引擎崩溃"这个桶都没进），而正常抛出时 `n_dock_failed=1`；
  同一场景两次运行给出**不同的漏斗台账**。在当前"报告未统一汇总降级标记"的背景下，这属于第 2 轮"降级台账"要治的同一类病。
- 第三处不一致：`slow-main` 场景下报告 `candidates=2` 而 checkpoint 已有 6 条 —— 报告会把**已经真实完成**的分子也漏掉。

缓解现状：`ScreeningFunnelReport.is_conserved`（`models.py:538`）会给出 `False`，日志里也打了
`"Screening report finalized: 6 final candidates. Conserved: False."`。
这是目前唯一的"事后可发现"手段，但它**不是异常**，也不阻止调用方继续。

---

## 6. 最小修复建议（产品侧；本轮不落盘）

### 6.1 推荐：把"硬中止"和"协作取消"拆成两个信号（3~6 行）

```python
# __init__ 里加: self._hard_abort: BaseException | None = None

def _do_single_dock(lig):
    if self._cancel_requested:
        raise KeyboardInterrupt("Docking cancelled")
    ...
    try:
        res = self._dock_single_safe(lig, box, params, sub_dir)
    except BaseException as exc:            # 1506
        self._hard_abort = exc              # ← 新增：留住原始异常
        self._cancel_requested = True
        raise
    ...

# 循环结束后（with 块之后、_finalize_report 之前）
if self._hard_abort is not None:            # ← 新增：与顺序无关地重放
    raise self._hard_abort
```

要点：`except` 分支先保存异常对象，**循环退出后再统一重放**，于是"主线程有没有先看到那个 Future"不再影响结果。
（也可以不加新字段，改成循环后用 `for f in future_to_lig: if f.done() and isinstance(f.exception(), BaseException): raise f.exception()` 收账——同样与顺序无关。）

### 6.2 一并建议（可选，防止再犯）

- `_run_single_tier` 在 `break` 之后、`_finalize_report` 之前，把"未处理的已崩溃 future"记账进 `n_dock_failed`，
  避免吞没/不吞没两条路径给出不同的漏斗台账。
- 给 `ScreeningFunnelReport` 增加显式的 `aborted: bool` + `abort_reason: str | None`，
  让"中止"成为一等公民，而不是只能靠 `is_conserved=False` 猜。

### 6.3 明确反对的"修法"

- ❌ 给用例加 sleep / 调等待（治不了机制，只是把窗口挪出观测范围；任务书也禁止）
- ❌ 放宽断言或删用例（等于把"崩溃必须抛出"这条契约删掉）
- ⚠️ 临时让 CI 绿只能 deselect 并挂缺陷单：`DOCKSTUDIO_GATES_DESELECT='tests/test_adversarial_m3_challenger2.py::test_mid_run_crash_simulation_at_7_and_19_zero_loss_resume'`，
  且必须注明"未修复，只是不阻塞门禁"。

### 6.4 修复后怎么验收（给下游）

```bash
# 缺陷是否还在：只要出现 NO_EXCEPTION 就说明还在（脚本退出码 0 = 仍可复现）
python scripts/qa/repro_flaky_resume.py --mode slow-main --repeat 10   # 修复前 10/10 吞
python scripts/qa/repro_flaky_resume.py --mode bare --repeat 20        # 修复前 ~35% 吞
# 用例本身
for i in $(seq 1 24); do python -m pytest tests/test_adversarial_m3_challenger2.py::test_mid_run_crash_simulation_at_7_and_19_zero_loss_resume -q; done
```
修复后的期望：`slow-main` 也 100% 抛出（因为重放与顺序无关），24 次隔离运行 24/24 通过。

---

## 7. 复现命令汇总

```bash
PY=/home/kui/miniconda3/envs/ai-gpu/bin/python
WT=/home/kui/wt/flaky          # 本分支的 worktree
LOCK=/home/kui/.dockstudio-heavy.lock

# 1) 失败率（24 次隔离）
flock $LOCK -c "cd $WT/server && for i in \$(seq 1 24); do \
  $PY -m pytest tests/test_adversarial_m3_challenger2.py::test_mid_run_crash_simulation_at_7_and_19_zero_loss_resume -q --tb=line; done"

# 2) 逐档对比（bare / slow-main / slow-worker / boundary）
flock $LOCK -c "cd $WT && $PY scripts/qa/repro_flaky_resume.py --all --repeat 10 -v"

# 3) 带日志看"静默返回"的那一次（会打印 Conserved: False）
flock $LOCK -c "cd $WT/server && $PY -m pytest \
  tests/test_adversarial_m3_challenger2.py::test_mid_run_crash_simulation_at_7_and_19_zero_loss_resume \
  -q -s -o log_cli=true --log-cli-level=INFO"
```

`scripts/qa/repro_flaky_resume.py` 的模式说明：`bare` 原样、`slow-main` 拖慢消费端、`slow-worker` 拖慢生产端、
`boundary` 打点（打印 `Future.result` / `_mark_dock_failed` / `_run_single_tier` / `run` 的事件序列）。
脚本只读产品代码，不写仓库任何文件；日志目录用 `/tmp`。

---

## 8. 边界与未做的事

- **没改产品代码**（本任务硬约束）：§6 的补丁只是文本建议，未落盘、未验证编译。
- **没改测试**：用例保持原样，本分支只加了一个取证脚本 + 本报告。
- **tiered 路径只做了 5 次对照**（`probe7`），结论"不吞"主要依据代码阅读；没有做 20+ 次的统计。
- **失败率的绝对数值随机器负载浮动**（同机队友的 gnina/redock 重活会改变线程时序）：62.5%（pytest 24 次）与 35%（脚本 20 次）都是"高频"，不是精确常数；能稳定复现的判据是 `slow-main` 模式下 100% 吞没。
- 未排查的相邻问题（另案）：三层路径崩溃后不排空队列（白跑完剩余分子）；`_mark_dock_failed` 在吞没路径漏记 `n_dock_failed`。
