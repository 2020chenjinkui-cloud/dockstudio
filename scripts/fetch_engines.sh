#!/usr/bin/env bash
# 下载各平台官方 AutoDock Vina 二进制到 resources/bin/（CI 与发版用）
set -euo pipefail

VERSION="${VINA_VERSION:-1.2.7}"
BASE="https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v${VERSION}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BIN="$ROOT/resources/bin"

fetch() {  # fetch <asset> <dest>
  local asset="$1" dest="$2"
  mkdir -p "$(dirname "$dest")"
  if [ -f "$dest" ]; then echo "exists: $dest"; return; fi
  echo "downloading $asset"
  curl -sfL -o "$dest" "$BASE/$asset"
  chmod +x "$dest"
}

fetch "vina_${VERSION}_win.exe"          "$BIN/windows/vina.exe"
fetch "vina_${VERSION}_linux_x86_64"     "$BIN/linux-x86_64/vina"
fetch "vina_${VERSION}_linux_aarch64"    "$BIN/linux-aarch64/vina"
fetch "vina_${VERSION}_mac_x86_64"       "$BIN/mac-x86_64/vina"
fetch "vina_${VERSION}_mac_aarch64"      "$BIN/mac-aarch64/vina"
echo "done. verify with: each binary --version"
