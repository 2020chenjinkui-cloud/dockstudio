#!/usr/bin/env python
"""崩溃续跑用例不稳定性的取证复现器（QA 工具，只读产品代码，不修改任何仓库文件）。

对应缺陷：tests/test_adversarial_m3_challenger2.py::
          test_mid_run_crash_simulation_at_7_and_19_zero_loss_resume
结论与证据链见 docs/qa/flaky-resume-2026-11-18.md。

四种模式（都用同一个 mock 引擎 / 同一份 30 配体库，只改变线程时序）：

  bare         原样跑（不加任何延迟）—— 观察天然失败率
  slow-main    主线程消费端在 _emit_progress 里 sleep 0.3s（放大「flag 已被 worker 置位」窗口）
  slow-worker  worker 在每次 _dock_single_safe 前 sleep 0.3s（把 worker 拖到消费端后面）
  boundary     边界打点：Future.result / _mark_dock_failed / _run_single_tier / run

用法（在仓库根目录或 server/ 下都可）：
    python scripts/qa/repro_flaky_resume.py --mode bare --repeat 20
    python scripts/qa/repro_flaky_resume.py --mode slow-main --repeat 3
    python scripts/qa/repro_flaky_resume.py --all --repeat 10

退出码：只要能复现出「NO_EXCEPTION」就返回 0（说明缺陷仍在）；全部正常抛出返回 1。
"""

from __future__ import annotations

import argparse
import collections
import sys
import tempfile
import threading
import time
import traceback
from concurrent.futures import Future
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SERVER_DIR = REPO_ROOT / "server"
sys.path.insert(0, str(SERVER_DIR))
sys.path.insert(0, str(SERVER_DIR / "tests"))

import test_adversarial_m3_challenger2 as T  # noqa: E402
from dockstudio.core.screen_runner import ScreeningRunner  # noqa: E402
from dockstudio.models import ScreeningConfig  # noqa: E402

CRASH_IDS = {"LIG_007", "LIG_019"}


def _scenario(tmp: Path):
    prep_lib, _ = T._build_test_library(tmp, n_ligands=30)
    rec = T._make_dummy_pdbqt(tmp / "receptor.pdbqt")
    engine = T.MockAdversarialCrashEngine(crash_at_ids=set(CRASH_IDS))
    runner = ScreeningRunner(
        config=ScreeningConfig(enable_tiered_funnel=False, cpu=2),
        out_dir=tmp / "run",
        receptor_pdbqt=rec,
        prepared_library=prep_lib,
        engine=engine,
        cpu_per_job=2,
    )
    return runner, engine


def _install_delays(mode: str, delay: float) -> None:
    orig_emit = ScreeningRunner._emit_progress
    orig_safe = ScreeningRunner._dock_single_safe

    def slow_emit(self, *a, **kw):
        if mode == "slow-main" and delay:
            time.sleep(delay)
        return orig_emit(self, *a, **kw)

    def slow_worker(self, lig, box, params, work_dir):
        if mode == "slow-worker" and delay:
            time.sleep(delay)
        return orig_safe(self, lig, box, params, work_dir)

    ScreeningRunner._emit_progress = slow_emit
    ScreeningRunner._dock_single_safe = slow_worker


def _install_boundary_hooks(events: list[str]) -> None:
    orig_result = Future.result
    orig_mark = ScreeningRunner._mark_dock_failed
    orig_st = ScreeningRunner._run_single_tier
    orig_run = ScreeningRunner.run

    def th() -> str:
        return threading.current_thread().name.replace("ThreadPoolExecutor-", "T")

    def spy_result(self, *a, **kw):
        try:
            return orig_result(self, *a, **kw)
        except BaseException as exc:
            events.append(f"Future.result RAISED {type(exc).__name__} in {th()}")
            raise

    def spy_mark(self, ligand_id):
        events.append(f"_mark_dock_failed({ligand_id}) in {th()}")
        return orig_mark(self, ligand_id)

    def spy_st(self, *a, **kw):
        try:
            r = orig_st(self, *a, **kw)
        except BaseException as exc:
            events.append(f"_run_single_tier RAISED {type(exc).__name__} in {th()}")
            raise
        events.append(f"_run_single_tier RETURNED in {th()}")
        return r

    def spy_run(self, *a, **kw):
        try:
            r = orig_run(self, *a, **kw)
        except BaseException as exc:
            events.append(f"run RAISED {type(exc).__name__} in {th()}")
            raise
        events.append(f"run RETURNED in {th()}")
        return r

    Future.result = spy_result
    ScreeningRunner._mark_dock_failed = spy_mark
    ScreeningRunner._run_single_tier = spy_st
    ScreeningRunner.run = spy_run


def run_mode(mode: str, repeat: int, delay: float, verbose: bool) -> int:
    _install_delays(mode, delay)
    events: list[str] = []
    if mode == "boundary":
        _install_boundary_hooks(events)

    tally: collections.Counter[str] = collections.Counter()
    for it in range(repeat):
        events.clear()
        tmp = Path(tempfile.mkdtemp(prefix=f"repro-{mode}-", dir="/tmp"))
        runner, engine = _scenario(tmp)
        try:
            report = runner.run(resume=False)
            outcome = "NO_EXCEPTION"
            tally[outcome] += 1
            detail = (
                f"candidates={len(report.candidates)} n_ready={report.n_ready} "
                f"n_docked_success={report.n_docked_success} n_dock_failed={report.n_dock_failed} "
                f"is_conserved={report.is_conserved}"
            )
        except BaseException as exc:  # noqa: BLE001 - 取证：任何抛出都算一次 PROPAGATE
            outcome = f"RAISED {type(exc).__name__}"
            tally["RAISED"] += 1
            detail = str(exc)
        ckpt = tmp / "run" / "checkpoint.jsonl"
        n_ck = len([ln for ln in ckpt.read_text().splitlines() if ln.strip()]) if ckpt.is_file() else -1
        print(
            f"  [{mode}] iter {it + 1:02d}: {outcome:22} dock_calls={len(engine.call_history):2} "
            f"checkpoint_lines={n_ck:2} | {detail}",
            flush=True,
        )
        if mode == "boundary" and verbose:
            for ev in events:
                print("        " + ev, flush=True)

    print(f"  [{mode}] 合计: {dict(tally)}（repeat={repeat}, delay={delay}s）", flush=True)
    return 0 if tally.get("NO_EXCEPTION") else 1


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mode", choices=["bare", "slow-main", "slow-worker", "boundary"], default="bare")
    ap.add_argument("--repeat", type=int, default=20)
    ap.add_argument("--delay", type=float, default=0.3)
    ap.add_argument("--all", action="store_true", help="依次跑 bare / slow-main / slow-worker / boundary")
    ap.add_argument("-v", "--verbose", action="store_true", help="boundary 模式打印事件序列")
    args = ap.parse_args()

    modes = ["bare", "slow-main", "slow-worker", "boundary"] if args.all else [args.mode]
    reproduced = False
    for mode in modes:
        print(f"\n=== mode={mode} ===", flush=True)
        rc = run_mode(mode, args.repeat, args.delay, args.verbose)
        reproduced = reproduced or rc == 0
    if reproduced:
        print("\n缺陷复现成功：至少一种模式下 run() 静默返回（NO_EXCEPTION），中止信号被吞。")
        return 0
    print("\n本次未复现：所有模式都正常抛出（说明缺陷可能已修复，或时序窗口未命中）。")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
