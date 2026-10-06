#!/usr/bin/env bash
# 构建后端冻结产物（本地调试用；CI 走 .github/workflows/ci.yml）
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

PY="python3"
if [[ "$OSTYPE" == msys* || "$OSTYPE" == cygwin* || -n "${WINDIR:-}" ]]; then
  PY="py -3.12"
fi

bash scripts/fetch_engines.sh
$PY -m pip install -e "server[dev]" pyinstaller
$PY -m PyInstaller scripts/backend.spec --noconfirm --clean
echo "output: $ROOT/dist/dockstudio-backend/"
