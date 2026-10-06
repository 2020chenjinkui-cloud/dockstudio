#!/usr/bin/env bash
# 把 PyInstaller 产物重命名为 Tauri sidecar 要求的 <name>-<target-triple>(.exe) 形式
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BIN_DIR="$ROOT/apps/desktop/src-tauri/binaries"
TRIPLE="$(rustc -vV | sed -n 's/host: //p')"

cd "$BIN_DIR"
if [[ "$TRIPLE" == *windows* ]]; then
  cp dockstudio-backend.exe "dockstudio-backend-${TRIPLE}.exe"
else
  cp dockstudio-backend "dockstudio-backend-${TRIPLE}"
fi
echo "staged sidecar: dockstudio-backend-${TRIPLE}"
