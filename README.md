# DockStudio

**[中文文档](README.zh-CN.md)**

Cross-platform molecular docking platform for **Windows / macOS / Linux** — friendly enough for first-time users, powerful enough for batch virtual screening veterans.

[![Release](https://img.shields.io/github/v/release/2020chenjinkui-cloud/dockstudio)](https://github.com/2020chenjinkui-cloud/dockstudio/releases)
[![CI](https://github.com/2020chenjinkui-cloud/dockstudio/actions/workflows/ci.yml/badge.svg)](https://github.com/2020chenjinkui-cloud/dockstudio/actions/workflows/ci.yml)
[![Platform](https://img.shields.io/badge/platform-Windows%20%7C%20macOS%20%7C%20Linux-blue)](https://github.com/2020chenjinkui-cloud/dockstudio/releases)

## Why DockStudio

Molecular docking usually means juggling CLI tools, file-format converters and scattered Python scripts. DockStudio wraps the whole pipeline — receptor/ligand preparation, pocket prediction, docking, interaction analysis and quality control — into one desktop app with **full provenance tracking**: every run can be audited, replayed and reproduced.

- **Simple mode**: a five-step wizard with a one-click built-in example (PDB: 1HSG) that takes you from zero to a finished docking result in minutes.
- **Advanced mode**: batch virtual screening with job queues, ranking, CSV export, flexible residues and an equivalent command line you can copy to reproduce any run outside the GUI.

## Features

### Docking engines
- **AutoDock Vina 1.2.7** bundled out of the box — nothing to install
- **smina** engine for Vinardo scoring & off-docking workflows
- **gnina / Uni-Dock** GPU engines as optional plug-ins (Linux + NVIDIA), with runtime path auto-detection
- Engine binaries ship per-platform in the installer; external engines register via settings

### Two working modes
- Five-step wizard: receptor → ligand → pocket → docking → analysis, each step validated with clear guidance
- Built-in **one-click example project**: ligand/receptor of HIV-1 protease (1HSG) pre-staged for a full end-to-end run
- Batch virtual screening: molecule library queues, pause/cancel, affinity ranking, one-click CSV export
- Flexible receptor residues support with plain-English validation errors
- **Equivalent CLI reproduction**: every GUI run shows the exact `vina`/`smina` command it executes

### Analysis & quality control
- 3D visualization (3Dmol.js) of receptor, pocket and poses
- **PLIP** interaction profiling (H-bonds, hydrophobic, π-stacking, salt bridges…)
- 2D interaction diagrams, screenshot export
- **PoseBusters** pose sanity checks
- Co-crystal ligand RMSD self-docking validation

### Provenance & reproducibility (v2.0)
- Full **event stream (JSONL)** for every step of a run
- One-click **`audit.zip`** audit bundle: inputs (with SHA256), configs, results, environment fingerprint (rdkit / meeko / openbabel versions), engine versions
- **`dockstudio audit` CLI** to inspect and re-verify audit bundles — agent/script friendly
- Pinned dependency upper bounds (e.g. `rdkit<2026.4`) for cross-environment reproducibility

### Getting the right input data
- In-app **input format guidance cards**: what a valid receptor PDB / ligand SDF / SMILES looks like
- Source wizards pointing to **RCSB PDB**, **PubChem** and **ZINC** with step-by-step download instructions
- SMILES quick lookup + built-in **demo molecule library** (10 real approved drugs)
- Bilingual UI (中文 / English)

## Download & install

Grab the latest installer from [Releases](https://github.com/2020chenjinkui-cloud/dockstudio/releases):

| Platform | File | Notes |
|---|---|---|
| Windows | `DockStudio_*_x64-setup.exe` (NSIS) or `*_x64_en-US.msi` | Per-user install, no admin needed; WebView2 required (preinstalled on Win11) |
| macOS (Apple Silicon) | `DockStudio_*_aarch64.dmg` | First launch: right-click → Open (unsigned) |
| Linux | `*.AppImage` / `*_amd64.deb` / `*.x86_64.rpm` | AppImage needs FUSE (install `libfuse2` on modern distros) |

GPU engines (gnina / Uni-Dock) currently support **Linux + NVIDIA**; the installer ships CPU engines for all platforms, GPU engines are detected from your `PATH` / engine settings.

## Quick start

1. Install & launch DockStudio.
2. Choose **Simple mode** → click **Load example project (1HSG)**.
3. Walk the five steps with defaults — you'll get a docked pose, affinity score, PLIP interactions and a PoseBusters report.
4. Open **Export → audit bundle** to see the provenance package, or copy the **equivalent command line** and re-run it in your terminal.

Running your own targets? Drop a receptor PDB (from RCSB) and a ligand (SDF from PubChem, or paste a SMILES) — the guidance cards in each step explain the accepted formats.

## Docker

```bash
docker compose up --build
# Backend API at http://127.0.0.1:8000 (health: /api/health)
```

The image ships the backend in **server mode** (API only, no UI) — consume it via the REST API / SDK, or point a desktop client at it. An optional GPU worker profile (`docker compose --profile gpu up`) runs a second instance with gnina / Uni-Dock enabled (requires the NVIDIA Container Toolkit).

## Development

```bash
# Backend (FastAPI + domain core + engine plugins, package name: dockstudio)
cd server
pip install -e ".[dev]"
uvicorn dockstudio.api.app:create_app --factory --reload

# Frontend (React 18 + Vite + TypeScript + 3Dmol.js)
cd apps/web
npm install
npm run dev

# Desktop shell (Tauri 2, Rust + sidecar backend)
cd apps/desktop && npm install
npm run tauri dev
```

CLI entry points after `pip install -e server`: `dockstudio` (incl. `dockstudio audit`) and `dockstudio-server`.

### Repository layout

```
server/      Python backend: FastAPI API + core domain + engine plugins
apps/web/    React + Vite frontend (3Dmol.js viewer)
apps/desktop/ Tauri 2 desktop shell (Rust)
resources/   Engine binaries shipped in installers (per-platform, not in git)
scripts/     Packaging, smoke-test, release scripts (fetch_engines.sh, pyinstaller specs…)
docker/      Dockerfile + docker-compose.yml
docs/        Dev docs: API contract, audit schema, SDK, deploy, AI-agent guide
```

## Testing & QA

- **123 backend tests** (pytest, incl. real-Vina end-to-end) + ruff lint
- **98 frontend tests**
- **Golden redocking regression** over 10 co-crystal targets (weekly self-hosted CI)
- Installers are smoke-tested for real: silent NSIS install on Windows, DEB/AppImage on physical Linux machines
- CI: GitHub Actions (GitHub-hosted for win/mac, self-hosted runner for Linux builds & GPU-adjacent jobs)

## Known limitations

- **Boltz-2 AI affinity re-scoring** is supported as an optional step but requires your own `boltz` installation and a GPU with ≥10 GB VRAM; verified hardware testing is still pending
- macOS packages are built & CI-tested but have not yet been validated on a physical Mac
- gnina / Uni-Dock GPU engines: Linux + NVIDIA only

## Documentation

- [API contract](docs/api-contract.md) · [Audit bundle schema](docs/audit-schema.md) · [SDK](docs/sdk.md) · [Deployment](docs/deploy.md) · [AI-agent guide](docs/ai.md)

## License

License not yet declared — all rights reserved by the author. If you plan to reuse code, please open an issue to discuss.
