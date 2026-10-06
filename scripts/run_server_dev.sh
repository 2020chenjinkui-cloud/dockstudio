#!/usr/bin/env bash
# Start the DockStudio dev server with repo-local data dir.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"

export DOCKSTUDIO_DATA_DIR="${DOCKSTUDIO_DATA_DIR:-$REPO_ROOT/data}"
mkdir -p "$DOCKSTUDIO_DATA_DIR"

cd "$REPO_ROOT/server"
exec py -3.12 -m uvicorn dockstudio.api.app:create_app --factory --reload \
    --host 127.0.0.1 --port 8000
