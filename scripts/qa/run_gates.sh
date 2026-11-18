#!/usr/bin/env bash
# DockStudio 2.0.5 — 一键质量门禁 / one-command quality gate
#
# 级别（逐级包含）:
#   --level l0   lint 门：ruff check + ruff format --check           （秒级）
#   --level l1   l0 + 后端快速子集 pytest（排除真机 e2e / golden）   （分钟级，默认）
#   --level l2   l0 + 全量后端 pytest（含真机 vina e2e）             （十几分钟）
#
# 用法:
#   scripts/qa/run_gates.sh                 # 默认 l1
#   scripts/qa/run_gates.sh --level l0      # 只跑 lint
#   scripts/qa/run_gates.sh --level l2      # 全量
#   scripts/qa/run_gates.sh --skip-lint     # 跳过 l0 两步（任何 level）
#   scripts/qa/run_gates.sh --help
#
# 任意一步失败 → 跑完全部步骤后以非零退出（=门禁红）。
#
# 环境变量（都有默认值，按需覆盖）:
#   DOCKSTUDIO_PY              python 解释器
#                              （默认 /home/kui/miniconda3/envs/ai-gpu/bin/python，不存在则回退 python3）
#   DOCKSTUDIO_RUFF            ruff 命令行（默认 "<DOCKSTUDIO_PY> -m ruff"）
#   DOCKSTUDIO_HEAVY_LOCK      重活锁文件（默认 /home/kui/.dockstudio-heavy.lock；l1/l2 自动 flock）
#   DOCKSTUDIO_GATES_NO_LOCK=1 不用 flock（不建议：会和同机其它重活抢 CPU/GPU）
#   DOCKSTUDIO_GATES_NO_ENGINE=1 不自动注入 DOCKSTUDIO_VINA / DOCKSTUDIO_GNINA
#   DOCKSTUDIO_GATES_LOG_DIR   每步日志目录（默认 /tmp/dockstudio-gates）
#   DOCKSTUDIO_GATES_DESELECT  额外 deselect 的 node id（空格分隔，可跳过已知 race 用例）
#   DOCKSTUDIO_GATE_PYTEST_ARGS 追加给 pytest 的额外参数（空格分隔）
#
# 事实基线（数字出处）：docs/qa/baseline-2026-11-18.md
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
SERVER_DIR="$REPO_ROOT/server"

LEVEL="l1"
SKIP_LINT=0

# l1 排除的“重活”文件：真机引擎 e2e、作业生命周期 e2e、模块级真机 fixture。
# 名单来自实测 durations（见 docs/qa/baseline-2026-11-18.md），不是猜测。
L1_IGNORE=(
  "tests/test_cli_agent.py"          # 真机 demo prep+dock，实测 115s + 68s
  "tests/test_job_lifecycle_e2e.py"  # 真机作业生命周期 e2e，实测 63s
  "tests/test_e2e_1hsg.py"           # 真机 1HSG redocking e2e，实测 63s
  "tests/test_audit.py"              # 真机 dock + 审计包，实测 62s
  "tests/test_interactions.py"       # 模块级真机对接 fixture，实测 setup 61s
  "tests/test_posecheck.py"          # 模块级真机姿态 fixture，实测 setup 61s
  "tests/test_diagram.py"            # 模块级真机对接 fixture，实测 setup 61s
  "tests/test_flex.py"               # 柔性对接 e2e，实测 36s
  "tests/test_htvs_funnel_stress.py" # 高通量漏斗压测，实测 26s + 8s
  "tests/test_batch.py"              # 批量真机 e2e
)

usage() {
  cat <<'USAGE'
用法: scripts/qa/run_gates.sh [--level l0|l1|l2] [--skip-lint] [--help]

  --level l0   ruff check server/dockstudio + ruff format --check server/dockstudio
  --level l1   l0 + 后端快速 pytest 子集（排除真机 e2e / golden 集），默认级别
  --level l2   l0 + 全量后端 pytest（含真机 vina e2e，走重活锁）
  --skip-lint  跳过 l0 的两步（任何级别）
  --help       显示本帮助

退出码: 0 = 全绿；1 = 有步骤失败；2 = 参数/环境错误。
事实基线见 docs/qa/baseline-2026-11-18.md。
USAGE
  exit 0
}

die() { printf 'run_gates: %s\n' "$*" >&2; exit 2; }

while [ $# -gt 0 ]; do
  case "$1" in
    --level)      [ $# -ge 2 ] || die "--level 需要一个参数 (l0|l1|l2)"; LEVEL="$2"; shift 2 ;;
    --level=*)    LEVEL="${1#*=}"; shift ;;
    --skip-lint)  SKIP_LINT=1; shift ;;
    -h|--help)    usage ;;
    *)            die "未知参数: $1（见 --help）" ;;
  esac
done

case "$LEVEL" in
  l0|0) LEVEL="l0" ;;
  l1|1) LEVEL="l1" ;;
  l2|2) LEVEL="l2" ;;
  *)    die "--level 只能是 l0|l1|l2，收到: $LEVEL" ;;
esac

PY="${DOCKSTUDIO_PY:-/home/kui/miniconda3/envs/ai-gpu/bin/python}"
[ -x "$PY" ] || PY="$(command -v python3 || command -v python || true)"
[ -n "$PY" ] || die "找不到 python 解释器（用 DOCKSTUDIO_PY 指定）"
if [ -n "${DOCKSTUDIO_RUFF:-}" ]; then
  read -r -a RUFF_CMD <<< "$DOCKSTUDIO_RUFF"
else
  RUFF_CMD=( "$PY" -m ruff )
fi
LOCK="${DOCKSTUDIO_HEAVY_LOCK:-/home/kui/.dockstudio-heavy.lock}"
LOG_DIR="${DOCKSTUDIO_GATES_LOG_DIR:-/tmp/dockstudio-gates}"
LOCK_ENABLED=1
if [ "${DOCKSTUDIO_GATES_NO_LOCK:-0}" = "1" ] || ! command -v flock >/dev/null 2>&1; then
  LOCK_ENABLED=0
fi

mkdir -p "$LOG_DIR" || die "无法创建日志目录: $LOG_DIR"
STAMP="$(date +%Y%m%d-%H%M%S)"
RESULTS=()
FAIL_COUNT=0

run_step() {  # run_step <gate-name> <command...>
  local name="$1"; shift
  local log="$LOG_DIR/$STAMP-$name.log"
  local start=$SECONDS rc dur
  printf '\n=== [%s] %s\n' "$name" "$(date -Is)"
  printf '    $ %s\n' "$*"
  "$@" 2>&1 | tee "$log"
  rc=${PIPESTATUS[0]}
  dur=$((SECONDS - start))
  if [ "$rc" -eq 0 ]; then
    printf -- '--- [%s] PASS  (%ss)\n' "$name" "$dur"
    RESULTS+=("PASS  $name  ${dur}s")
  else
    printf -- '--- [%s] FAIL  rc=%s  (%ss)  日志: %s\n' "$name" "$rc" "$dur" "$log"
    RESULTS+=("FAIL  $name  ${dur}s  rc=$rc  log=$log")
    FAIL_COUNT=$((FAIL_COUNT + 1))
  fi
  return 0   # 不中断：跑完所有步骤再统一汇总
}

# 引擎注入：CI 靠 scripts/fetch_engines.sh 把 vina 放进 resources/bin/；
# 干净 worktree 没有该目录（.gitignore 里是 resources/bin/），真机用例会失效。
inject_engines() {
  if [ "${DOCKSTUDIO_GATES_NO_ENGINE:-0}" = "1" ]; then
    echo "引擎注入:    已禁用 (DOCKSTUDIO_GATES_NO_ENGINE=1)"
    return 0
  fi
  local d
  for d in linux-x86_64 linux-aarch64 mac-aarch64 mac-x86_64 windows; do
    if [ -z "${DOCKSTUDIO_VINA:-}" ] && [ -f "$REPO_ROOT/resources/bin/$d/vina" ]; then
      export DOCKSTUDIO_VINA="$REPO_ROOT/resources/bin/$d/vina"
    fi
    if [ -z "${DOCKSTUDIO_GNINA:-}" ] && [ -f "$REPO_ROOT/resources/bin/$d/gnina" ]; then
      export DOCKSTUDIO_GNINA="$REPO_ROOT/resources/bin/$d/gnina"
    fi
  done
  echo "引擎注入:    DOCKSTUDIO_VINA=${DOCKSTUDIO_VINA:-<未设置>}"
  echo "引擎注入:    DOCKSTUDIO_GNINA=${DOCKSTUDIO_GNINA:-<未设置>}"
}

fingerprint() {
  echo "仓库根:      $REPO_ROOT"
  echo "git HEAD:    $(git -C "$REPO_ROOT" rev-parse HEAD 2>/dev/null || echo '<非 git 仓库>')"
  echo "git 分支:    $(git -C "$REPO_ROOT" rev-parse --abbrev-ref HEAD 2>/dev/null || echo '-')"
  echo "python:      $PY -> $("$PY" -V 2>&1)"
  echo "ruff:        $("${RUFF_CMD[@]}" --version 2>&1)"
  echo "pytest:      $("$PY" -m pytest --version 2>&1 | head -1)"
  if command -v nvidia-smi >/dev/null 2>&1; then
    echo "GPU:         $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>/dev/null | head -1)"
  fi
  echo "uname:       $(uname -sr)"
  echo "重活锁:      $([ "$LOCK_ENABLED" = "1" ] && echo "$LOCK（flock）" || echo '关闭')"
  echo "门禁级别:    $LEVEL$([ "$SKIP_LINT" = "1" ] && echo '（--skip-lint）')"
  inject_engines
}

gate_lint() {
  run_step "l0-ruff-check"  "${RUFF_CMD[@]}" check server/dockstudio
  run_step "l0-ruff-format" "${RUFF_CMD[@]}" format --check server/dockstudio
}

run_pytest() {  # run_pytest <gate-name> <pytest args...>（在 server/ 下跑，l1/l2 走重活锁）
  local name="$1"; shift
  local pre=()
  [ "$LOCK_ENABLED" = "1" ] && pre=( flock "$LOCK" )
  run_step "$name" bash -c 'cd "$1" || exit 1; shift; exec "$@"' \
    _ "$SERVER_DIR" "${pre[@]}" "$PY" -m pytest "$@"
}

gate_pytest_l1() {
  local args=( tests -q -m "not golden" )
  local f d
  for f in "${L1_IGNORE[@]}"; do args+=( "--ignore=$f" ); done
  for d in ${DOCKSTUDIO_GATES_DESELECT:-}; do args+=( "--deselect=$d" ); done
  read -r -a _extra <<< "${DOCKSTUDIO_GATE_PYTEST_ARGS:-}"
  [ "${#_extra[@]}" -gt 0 ] && args+=( "${_extra[@]}" )
  echo "l1 排除的重活文件（实测 durations 见 docs/qa/baseline-2026-11-18.md）："
  for f in "${L1_IGNORE[@]}"; do echo "  - $f"; done
  run_pytest "l1-pytest-fast" "${args[@]}"
}

gate_pytest_l2() {
  local args=( tests -q )
  local d
  for d in ${DOCKSTUDIO_GATES_DESELECT:-}; do args+=( "--deselect=$d" ); done
  read -r -a _extra <<< "${DOCKSTUDIO_GATE_PYTEST_ARGS:-}"
  [ "${#_extra[@]}" -gt 0 ] && args+=( "${_extra[@]}" )
  run_pytest "l2-pytest-full" "${args[@]}"
}

cd "$REPO_ROOT" || die "无法进入仓库根: $REPO_ROOT"

echo "############ DockStudio QA gates ############"
fingerprint

[ "$SKIP_LINT" = "1" ] || gate_lint
case "$LEVEL" in
  l0) : ;;
  l1) gate_pytest_l1 ;;
  l2) gate_pytest_l2 ;;
esac

echo
echo "############ 汇总 (level=$LEVEL) ############"
for r in "${RESULTS[@]}"; do echo "  $r"; done
echo "日志目录: $LOG_DIR"
if [ "$FAIL_COUNT" -eq 0 ]; then
  echo "门禁结果: 绿（0 个失败步骤）"
  exit 0
else
  echo "门禁结果: 红（$FAIL_COUNT 个失败步骤）"
  exit 1
fi
