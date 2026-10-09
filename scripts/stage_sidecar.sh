#!/usr/bin/env bash
# 把 onefile 冻结后端（dist/dockstudio-backend(.exe)）重命名为
# Tauri sidecar 要求的 binaries/<name>-<target-triple>(.exe) 形式
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BIN_DIR="$ROOT/apps/desktop/src-tauri/binaries"
mkdir -p "$BIN_DIR"
TRIPLE="$(rustc -vV | sed -n 's/host: //p')"

if [[ "$TRIPLE" == *windows* ]]; then
  SRC="$ROOT/dist/dockstudio-backend.exe"
  DST="$BIN_DIR/dockstudio-backend-${TRIPLE}.exe"
else
  SRC="$ROOT/dist/dockstudio-backend"
  DST="$BIN_DIR/dockstudio-backend-${TRIPLE}"
fi
[[ -f "$SRC" ]] || { echo "onefile backend missing: $SRC （先跑 backend_onefile.spec）"; exit 1; }
cp "$SRC" "$DST"
ls -la "$DST"
echo "staged sidecar: dockstudio-backend-${TRIPLE}"
