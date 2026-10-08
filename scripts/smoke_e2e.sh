#!/usr/bin/env bash
# 端到端冒烟：起后端 → 跑通 1HSG 真实对接 → 断言 top affinity 区间 → 停服务。
# 用法:
#   bash scripts/smoke_e2e.sh [后端可执行路径]   本地：dev 模式（uvicorn）或冻结产物
#   bash scripts/smoke_e2e.sh --base-url http://host:port   对已运行的服务做冒烟（远端/容器）
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DATA_DIR="$(mktemp -d)"
FIXTURES="$ROOT/server/tests/fixtures"
export DOCKSTUDIO_DATA_DIR="$DATA_DIR"

PY="python3"
if [[ -n "${WINDIR:-}" ]]; then PY="py -3.12"; fi

BASE_URL=""
BACKEND=""
if [[ "${1:-}" == "--base-url" ]]; then
  BASE_URL="${2:?missing url}"
else
  BACKEND="${1:-}"
fi

# 冻结后端 main() 固定绑 8000；dev 模式自己起 uvicorn 才用随机端口
if [[ -n "${PORT:-}" ]]; then
  :
elif [[ -n "$BACKEND" ]]; then
  PORT=8000
else
  PORT=$((20000 + RANDOM % 20000))
fi

PID=""
LISTEN_PID=""
cleanup() {
  # py.exe 启动器被杀后 python.exe 子进程会变孤儿；按监听端口补杀
  if [[ -n "$LISTEN_PID" ]]; then taskkill //F //T //PID "$LISTEN_PID" 2>/dev/null || true; fi
  [[ -n "$PID" ]] && kill "$PID" 2>/dev/null || true
}
trap cleanup EXIT

if [[ -z "$BASE_URL" ]]; then
  if [[ -n "$BACKEND" ]]; then
    "$BACKEND" & PID=$!
  else
    cd "$ROOT/server"
    DOCKSTUDIO_DATA_DIR="$DATA_DIR" $PY -m uvicorn dockstudio.api.app:create_app --factory --port "$PORT" &
    PID=$!
    cd "$ROOT"
  fi
fi

BASE="${BASE_URL:-http://127.0.0.1:$PORT}"
for i in $(seq 1 60); do
  curl -sf "$BASE/api/health" > /dev/null 2>&1 && break
  sleep 1
  if [[ "$i" == 60 ]]; then echo "server did not start"; exit 1; fi
done
LISTEN_PID="$(netstat -ano 2>/dev/null | grep ":$PORT" | grep LISTEN | awk '{print $NF}' | head -1 || true)"

echo "[smoke] engines:"
curl -sf "$BASE/api/engines" | head -c 200; echo

curl -sf -F "file=@$FIXTURES/1HSG.pdb" "$BASE/api/prepare/receptor" > "$DATA_DIR/rec.json"
curl -sf -F "smiles=$(cat "$FIXTURES/indinavir_smiles.txt")" -F "name=indinavir" "$BASE/api/prepare/ligand" > "$DATA_DIR/lig.json"

RESULT=$($PY - <<EOF
import json, os, subprocess, sys
data_dir = os.environ["DOCKSTUDIO_DATA_DIR"]
rec = json.load(open(os.path.join(data_dir, "rec.json"), encoding="utf-8"))
lig = json.load(open(os.path.join(data_dir, "lig.json"), encoding="utf-8"))
center = rec["cocrystal_ligands"][0]["center"]
body = {
    "engine": "vina",
    "receptor_id": rec["id"],
    "ligand_id": lig["id"],
    "box": {"center": center, "size": [22.5, 22.5, 22.5]},
    "params": {"exhaustiveness": 4, "num_modes": 5, "energy_range": 3.0,
               "cpu": 8, "seed": None, "scoring": "vina"},
}
out = subprocess.run(
    ["curl", "-sf", "-X", "POST", "$BASE/api/dock",
     "-H", "Content-Type: application/json", "-d", json.dumps(body)],
    capture_output=True, text=True)
if out.returncode != 0:
    sys.exit("dock request failed: " + out.stderr[:500])
print(out.stdout)
EOF
)

echo "$RESULT" | $PY -c "
import json, sys
r = json.load(sys.stdin)
poses = r['poses']
assert len(poses) >= 3, f'too few poses: {len(poses)}'
top = poses[0]['affinity']
assert -14.0 <= top <= -6.0, f'top affinity out of range: {top}'
print(f'[smoke] PASS  poses={len(poses)}  top_affinity={top} kcal/mol')
"
