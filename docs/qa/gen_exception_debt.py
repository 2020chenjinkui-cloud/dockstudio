#!/usr/bin/env python3
"""Generate docs/qa/exception-debt.md from the live tree.

Registers every BLE001/S110/S112 site of server/dockstudio with:
  file:line, module, enclosing function, handler shape, whether the module sits on a
  path that shapes scientific conclusions, and the round-2 priority.

Run from the repo root:  python3 /tmp/gen_exception_debt.py
"""

from __future__ import annotations

import ast
import collections
import os
import pathlib
import subprocess
import sys

REPO = pathlib.Path(__file__).resolve().parents[2]  # <repo>/docs/qa/gen_exception_debt.py
SERVER = REPO / "server"
PKG = SERVER / "dockstudio"
RUFF = os.environ.get("RUFF", "ruff")
BASELINE = os.environ.get("BASELINE_REF", "HEAD")
OUT = REPO / "docs/qa/exception-debt.md"

# ---------------------------------------------------------------- diagnostics
proc = subprocess.run(
    [RUFF, "check", "dockstudio", "--select", "BLE001,S110,S112", "--output-format", "concise"],
    cwd=SERVER,
    capture_output=True,
    text=True,
)
lines = [ln for ln in proc.stdout.splitlines() if ln.startswith("dockstudio/")]
diags: list[tuple[str, int, str]] = []
for ln in lines:
    head, code = ln.split(": ", 1)
    path, lineno, _col = head.rsplit(":", 2)
    diags.append((path, int(lineno), code.split()[0]))
if not diags:
    sys.exit(f"no diagnostics parsed; stdout={proc.stdout[:400]!r} stderr={proc.stderr[:400]!r}")

# ---------------------------------------------------------------- AST scanning
SCIENCE_DIRS = {"prep", "analysis", "core"}  # 直接决定数值结论 / 筛选结果
INDIRECT_DIRS = {"engines"}  # 引擎失败降级直接影响构象产出


class Handler:
    __slots__ = ("path", "node", "qualname", "line", "end")

    def __init__(self, path, node, qualname):
        self.path = path
        self.node = node
        self.qualname = qualname
        self.line = node.lineno
        self.end = node.end_lineno


def collect_handlers(path: pathlib.Path, rel: str) -> list[Handler]:
    return collect_handlers_from_source(path.read_text(encoding="utf-8"), rel)


def collect_handlers_from_source(src: str, rel: str) -> list[Handler]:
    tree = ast.parse(src)
    found: list[Handler] = []

    def walk(node, stack):
        for child in ast.iter_child_nodes(node):
            name = None
            if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef)):
                name = child.name
            elif isinstance(child, ast.ClassDef):
                name = child.name
            if name is not None:
                walk(child, stack + [name])
            else:
                if isinstance(child, ast.ExceptHandler):
                    found.append(Handler(rel, child, ".".join(stack) or "<module>"))
                walk(child, stack)

    walk(tree, [])
    return found


def body_kind(handler: ast.ExceptHandler) -> tuple[str, str]:
    """(kind, chinese label) for the except body."""
    body = handler.body
    src = handler
    if len(body) == 1 and isinstance(body[0], ast.Pass):
        return "pass", "except: pass（完全吞掉，无日志无痕迹）"
    if len(body) == 1 and isinstance(body[0], ast.Continue):
        return "continue", "except: continue（跳过当前项继续）"
    if any(isinstance(n, ast.Raise) for n in ast.walk(src)) and not any(
        isinstance(n, ast.Return) for n in body
    ):
        first = body[0]
        if isinstance(first, ast.Raise):
            return "raise", "转译为具体异常（不静默）"
    text = ast.unparse(handler)
    has_log = any(
        isinstance(n, ast.Call)
        and isinstance(n.func, ast.Attribute)
        and n.func.attr in {"warning", "info", "debug", "error", "exception", "critical"}
        for n in ast.walk(src)
    )
    has_warn = any(
        isinstance(n, ast.Call)
        and isinstance(n.func, ast.Attribute)
        and n.func.attr in {"append", "add", "extend"}
        for n in ast.walk(src)
    )
    if has_log:
        return "logged", "有日志记录（降级可见）"
    if has_warn:
        return "recorded", "记入 warnings/结果列表（降级可见）"
    if any(isinstance(n, ast.Return) for n in body):
        return "return", "返回降级默认值（无日志）"
    if "pass" in text and len(body) > 1:
        return "silent", "无可视痕迹（吞掉异常）"
    return "other", "吞掉异常后继续执行（无日志）"


def priority(rel: str, kind: str) -> tuple[str, str]:
    parts = rel.split("/")
    top = parts[1] if len(parts) > 1 else ""
    in_science = top in SCIENCE_DIRS
    if in_science and kind in {"pass", "continue", "return", "silent", "other"}:
        return "P0", "是"
    if in_science:
        return "P1", "是"
    if top in INDIRECT_DIRS:
        return "P1", "是"
    if kind in {"pass", "continue"}:
        return "P1", "否"
    return "P2", "否"


handlers: list[Handler] = []
for py in sorted(PKG.rglob("*.py")):
    rel = py.relative_to(SERVER).as_posix()
    handlers.extend(collect_handlers(py, rel))

# map diagnostics onto handlers
by_loc = collections.defaultdict(set)
unmapped = []
for rel, lineno, code in diags:
    hit = [h for h in handlers if h.path == rel and h.line <= lineno <= h.end]
    if not hit:
        unmapped.append((rel, lineno, code))
        continue
    if len(hit) > 1:  # nested handlers: take the innermost
        hit.sort(key=lambda h: h.end - h.line)
    by_loc[(rel, hit[0].line)].add(code)
if unmapped:
    sys.exit(f"unmapped diagnostics ({len(unmapped)}): {unmapped[:5]}")

# which sites were masked before this round (noqa / per-file-ignore)?
masked: dict[tuple[str, int], str] = {}
for rel in sorted({h.path for h in handlers}):
    cur = [h for h in handlers if h.path == rel]
    head_src = subprocess.run(
        ["git", "show", f"{BASELINE}:server/{rel}"], cwd=REPO, capture_output=True, text=True
    ).stdout
    if not head_src:
        continue
    head_lines = head_src.splitlines()
    head_handlers = collect_handlers_from_source(head_src, rel)
    assert len(head_handlers) == len(cur), (rel, len(head_handlers), len(cur))
    for h, old in zip(cur, head_handlers):
        old_line = head_lines[old.line - 1] if old.line - 1 < len(head_lines) else ""
        if "noqa: BLE001" in old_line:
            masked[(rel, h.line)] = "原 `# noqa: BLE001`（本轮改写为 `# 降级兜底` 注释）"
        elif "# noqa" in old_line:
            masked[(rel, h.line)] = "原 `# noqa`（本轮已移除）"
        if rel.endswith("dockstudio/__init__.py") and by_loc.get((rel, h.line), set()) & {"S110", "S112"}:
            prev = masked.get((rel, h.line))
            masked[(rel, h.line)] = (
                (prev + " + " if prev else "") + "原被 `per-file-ignores` 豁免 S110/S112（本轮移除）"
            )

# ---------------------------------------------------------------- roll-up
rows = []
for h in handlers:
    codes = by_loc.get((h.path, h.line))
    if not codes:
        continue
    kind, label = body_kind(h.node)
    prio, science = priority(h.path, kind)
    one_line = " ".join(ast.unparse(h.node).split())
    if len(one_line) > 76:
        one_line = one_line[:73] + "..."
    rows.append(
        {
            "path": h.path,
            "line": h.line,
            "module": "/".join(h.path.split("/")[:-1]),
            "qualname": h.qualname,
            "codes": "+".join(sorted(codes)),
            "handler": one_line.replace("|", "\\|"),
            "kind": kind,
            "kind_label": label,
            "priority": prio,
            "science": science,
            "masked": masked.get((h.path, h.line), ""),
        }
    )

rows.sort(key=lambda r: (r["priority"], r["path"], r["line"]))
n_codes = collections.Counter()
for rel, _ln, code in diags:
    n_codes[code] += 1
n_sites = len(rows)
by_prio = collections.Counter(r["priority"] for r in rows)
by_kind = collections.Counter(r["kind"] for r in rows)
by_module = collections.Counter(r["module"] for r in rows)
science_rows = [r for r in rows if r["science"] == "是"]

# ---------------------------------------------------------------- emit
doc = []
w = doc.append


def table(rs, with_prio=True):
    w("| # | 位置 | 规则 | 所在函数 | 处理器形态 | 降级可见性 | 科学路径 | 第 2 轮优先级 |")
    w("|---|---|---|---|---|---|---|---|")
    for i, r in enumerate(rs, 1):
        w(
            f"| {i} | `{r['path'].removeprefix('dockstudio/')}:{r['line']}` | {r['codes']} | "
            f"`{r['qualname']}` | `{r['handler']}` | {r['kind_label']} | {r['science']} | "
            f"{r['priority']} |"
        )
    w("")


w("# 异常捕获债务登记（BLE001 / S110 / S112）")
w("")
w("> 本文档是 `ruff` 异常捕获类规则在 `server/dockstudio` 的**全量台账**。")
w("> 本轮（lint/format 门禁修复）**不修**这些代码，只把它们从静态检查里显式摘出来")
w("> 并逐条登记，作为第 2 轮“统一降级台账”的输入清单。")
w("")
w("| 口径 | 值 |")
w("|---|---|")
w("| 生成基线 | `main` @ c0f33cc（DockStudio 2.0.5） |")
w("| 工具 | ruff 0.16.10（`server/pyproject.toml` 已精确锁定） |")
w("| 规则分布 | BLE001 %d / S110 %d / S112 %d 条诊断 |" % (n_codes["BLE001"], n_codes["S110"], n_codes["S112"]))
w("| 去重后的 except 处理点 | **%d 处**（同一处 `except Exception: pass` 会同时命中 BLE001 与 S110/S112） |" % n_sites)
w("| 其中科学结论路径 | %d 处（P0 %d / P1 %d） |" % (
    len(science_rows),
    sum(1 for r in science_rows if r["priority"] == "P0"),
    sum(1 for r in science_rows if r["priority"] == "P1"),
))
w("| 门禁处置 | `BLE001`/`S110`/`S112` **不在** `[tool.ruff.lint].select` 里，由 `ignore` 显式注明理由；**零 `# noqa` 掩盖** |")
w("")

w("## 1. 为什么本轮不修，而是登记")
w("")
w("这三条规则命中 208 条诊断、%d 个处理点，但它们**不是**代码风格噪音：" % n_sites)
w("")
w("- 绝大多数是**刻意的降级兜底**：GPU / 对接引擎 / 可选插件（openbabel、plip、posebusters）")
w("  不可用时回退到 CPU 或内置几何算法，这是产品行为的一部分，删掉兜底会让整个流程在")
w("  缺依赖的机器上直接崩。")
w("- 机械地把 `except Exception` 收窄成具体异常类型**会改变降级行为**，而某个调用点到底")
w("  该不该降级、降级后要不要在报告里留痕，必须逐个按科学语义判断。")
w("- 反过来，一味 `# noqa` 或把整条规则关掉，正是文档里被反复复核抓住的")
w("  **“静默科学降级”** 的代码根源：降级发生了，但结果文件、报告和审计里看不到。")
w("")
w("所以本轮的处置是：**规则不进屋檐（不 select），代码不掩盖（不留 noqa），债务进台账**。")
w("")
w("## 2. 口径说明：175 条 vs 208 条")
w("")
w("修复前的 `ruff` 只报出 175 条（BLE001 133 / S110 37 / S112 5），比代码里真实存在的")
w("208 条少 33 条，差额全部来自**掩盖**：")
w("")
w("- 32 处 `except Exception:  # noqa: BLE001`（t1 基线里那句“168 处 except Exception 仅")
w("  32 处标注”就是这批）—— 本轮已把 noqa 改写成普通注释 `# 降级兜底`，于是它们重新")
w("  出现在统计里，因此 BLE001 从 133 涨到 165。")
w("- 1 处 `dockstudio/__init__.py` 里的 `try/except/pass|continue` 被 `per-file-ignores`")
w("  的 `S110/S112` 豁免（本轮一并移除），因此 S112 从 5 涨到 6。")
w("")
w("结论：**175 是“被压制后的报数”，208 才是债务真值**。本文档登记 208 条 / %d 个处理点，" % n_sites)
w("并在第 4 节标出这 32 个处理点（对应 33 条诊断）解除掩盖的位置，供第 2 轮核对。")
w("")
w("## 3. 统计")
w("")
w("### 3.1 按规则")
w("")
w("| 规则 | 含义 | 诊断数 |")
w("|---|---|---|")
w("| BLE001 | 盲捕 `Exception`/`BaseException` | %d |" % n_codes["BLE001"])
w("| S110 | `try-except-pass`：捕获后直接 `pass` | %d |" % n_codes["S110"])
w("| S112 | `try-except-continue`：捕获后 `continue` | %d |" % n_codes["S112"])
w("")
w("### 3.2 按处理点所在模块（去重后）")
w("")
w("| 模块 | 处理点 | 其中科学路径 |")
w("|---|---|---|")
for mod, cnt in sorted(by_module.items(), key=lambda kv: (-kv[1], kv[0])):
    sci = sum(1 for r in rows if r["module"] == mod and r["science"] == "是")
    w(f"| `{mod}` | {cnt} | {sci} |")
w("")
w("### 3.3 按降级可见性")
w("")
w("| 处理点形态 | 处理点数 |")
w("|---|---|")
KIND_LABEL = {
    "pass": "except: pass（完全吞掉，无日志无痕迹）",
    "continue": "except: continue（跳过当前项继续）",
    "return": "返回降级默认值（无日志）",
    "silent": "无可视痕迹（吞掉异常）",
    "logged": "有日志记录（降级可见）",
    "recorded": "记入 warnings/结果列表（降级可见）",
    "raise": "转译为具体异常（不静默）",
    "other": "吞掉异常后继续执行（无日志）",
}
for kind, cnt in sorted(by_kind.items(), key=lambda kv: (-kv[1], kv[0])):
    w(f"| {KIND_LABEL.get(kind, kind)} | {cnt} |")
w("")
w("> 第 2 轮最该先动的是“完全吞掉（pass）”“跳过当前项（continue）”“吞掉后继续执行”")
w("> “返回降级默认值”这几类：它们让“降级发生过”这件事在任何输出里都不可见，")
w("> 正是文档里被反复抓住的“静默科学降级”的定义。")
w("")
w("## 4. 逐条登记（按第 2 轮优先级排序）")
w("")
w("优先级定义：")
w("")
w("- **P0**：位于科学结论路径（`prep/`、`analysis/`、`core/`）**且**降级不可见 —— 会直接、")
w("  无声地影响 RMSD / 打分 / 筛选淘汰等数值结论。第 2 轮必须优先纳入降级台账。")
w("- **P1**：位于科学路径但降级可记录，或位于 `engines/`（引擎失败降级直接影响构象产出），")
w("  或虽然不在科学路径但属 `pass`/`continue` 型（用户请求被静默丢弃）。")
w("- **P2**：CLI / API / 展示层，降级不影响已产出的数值。")
w("")
PRIO_TITLE = {"P0": "P0 — 科学路径 + 降级不可见", "P1": "P1 — 科学路径已记录 / 引擎层 / 静默丢弃", "P2": "P2 — 交互与展示层"}
for prio in ("P0", "P1", "P2"):
    rs = [r for r in rows if r["priority"] == prio]
    w(f"### 4.{('P0','P1','P2').index(prio) + 1} {PRIO_TITLE[prio]}（{len(rs)} 处）")
    w("")
    if rs:
        table(rs)
    else:
        w("（无）")
        w("")
w("### 4.4 本轮解除掩盖的位置（32 个处理点 / 33 条诊断，第 2 轮核对用）")
w("")
w("| 位置 | 规则 | 解除的掩盖形式 |")
w("|---|---|---|")
n_masked = 0
for r in sorted(rows, key=lambda r: (r["path"], r["line"])):
    if r["masked"]:
        n_masked += 1
        w(f"| `{r['path'].removeprefix('dockstudio/')}:{r['line']}` | {r['codes']} | {r['masked']} |")
if n_masked == 0:
    w("| — | — | — |")
w("")
w(f"（共 {n_masked} 处。）")
w("")
w("## 5. 第 2 轮建议：统一降级台账")
w("")
w("本轮只做了“登记”，第 2 轮要把台账**做进代码和报告**，建议的落地形态：")
w("")
w("1. **降级事件结构体**：`DowngradeEvent(code, where, severity, detail, affects_result)`，")
w("   由统一的 `record_downgrade()` 写入。")
w("2. **标记常量集中管理**：现有散落的 `HIGH_STRAIN_UNCONVERGED_FALLBACK`、")
w("   `POSECHECK_*`、`*_FALLBACK` 等字符串标记收敛到一处枚举（本轮扫描已经能看到")
w("   `analysis/scoring_rigor.py`、`analysis/posecheck.py`、`core/screen_runner.py`")
w("   里同一语义存在多种写法）。")
w("3. **P0 处理点必须带 `affects_result=True`**：报表 / 交付包 / 审计 events 里汇总输出，")
w("   使“本次运行发生了哪些降级、影响了哪些分子”成为一等输出，而不是藏在日志行里。")
w("4. **门禁侧留一道闸**：当 `affects_result` 的降级比例超过阈值时把结果标记为不完全复现，")
w("   而不是照常给出 RMSD/富集数字。")
w("5. **回归**：第 2 轮改完后，把 `BLE001/S110/S112` 重新加入 `select` 并让剩余处理点")
w("   **逐点显式**（收窄异常类型或调用 `record_downgrade()`），而不是整条忽略。")
w("")
w("代码里已留锚点：原先 32 处 `# noqa: BLE001` 改写为 `# 降级兜底`，可直接定位：")
w("")
w("```bash")
w("grep -rn 降级兜底 server/dockstudio")
w("```")
w("")
w("## 6. 复现与维护")
w("")
w("```bash")
w("# 全量清单（与本文档同源）")
w("cd server && ruff check dockstudio --select BLE001,S110,S112 --output-format concise")
w("")
w("# 规则分布")
w("cd server && ruff check dockstudio --select BLE001,S110,S112 --statistics")
w("")
w("# 重新生成本文档（脚本随本轮提交在 docs/qa/）")
w("python3 docs/qa/gen_exception_debt.py")
w("```")
w("")
w("> 说明：生成脚本与本台账同目录（`docs/qa/gen_exception_debt.py`），不依赖第 1 轮")
w("> `scripts/qa/` 基线产物，避免跨分支依赖。")
w("")
w("变更记录：")
w("")
w("| 日期 | 变更 |")
w("|---|---|")
w("| 2026-11-18 | 首次生成：208 条诊断 / %d 个处理点；32 个处理点解除掩盖 |" % n_sites)
w("")

OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text("\n".join(doc) + "\n", encoding="utf-8")
print(f"wrote {OUT} ({len(doc)} lines, {n_sites} sites, diags={len(diags)})")
print("codes:", dict(n_codes), "priority:", dict(by_prio), "masked:", n_masked)
