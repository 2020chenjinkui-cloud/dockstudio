# DockStudio

跨平台（Windows / macOS / Linux）分子对接平台。详见：

- `../分子对接平台规划设计书.md`（上级目录）
- `../开发任务分解.md`（上级目录）

## 目录

- `server/` — Python 后端：FastAPI + 核心域 + 引擎插件（包名 `dockstudio`）
- `apps/web/` — React + Vite 前端（Mol* 3D 查看器）
- `apps/desktop/` — Tauri 桌面壳（需要 Rust 工具链，M0 只留配置骨架）
- `resources/` — 随包分发的引擎二进制（按平台分目录，不入 git）
- `scripts/` — 打包、冒烟、发版脚本
- `docs/` — API 契约等开发文档

## 开发模式快速开始

```bash
# 后端
cd server
pip install -e ".[dev]"
uvicorn dockstudio.api.app:create_app --factory --reload

# 前端（另开终端）
cd apps/web
npm install
npm run dev
```
