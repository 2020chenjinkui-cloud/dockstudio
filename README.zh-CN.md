# DockStudio

**[English](README.md)**

跨平台（**Windows / macOS / Linux**）分子对接平台 —— 新手十分钟出结果，老手批量虚拟筛选也够硬核。

[![Release](https://img.shields.io/github/v/release/2020chenjinkui-cloud/dockstudio)](https://github.com/2020chenjinkui-cloud/dockstudio/releases)
[![CI](https://github.com/2020chenjinkui-cloud/dockstudio/actions/workflows/ci.yml/badge.svg)](https://github.com/2020chenjinkui-cloud/dockstudio/actions/workflows/ci.yml)
[![Platform](https://img.shields.io/badge/平台-Windows%20%7C%20macOS%20%7C%20Linux-blue)](https://github.com/2020chenjinkui-cloud/dockstudio/releases)

## 为什么做 DockStudio

传统分子对接要在命令行工具、格式转换器、散装 Python 脚本之间来回倒腾。DockStudio 把完整流水线 —— 受体/配体准备、口袋预测、对接、相互作用分析、质检 —— 收敛进一个桌面应用，并且**全程可溯源**：每一次运行都可审计、可回放、可复现。

- **简易模式**：五步向导 + 内置一键示例项目（PDB: 1HSG），从零到完整对接结果只要几分钟
- **高级模式**：批量虚拟筛选，任务队列、排行、CSV 导出、柔性残基，并给出等价命令行供脱离 GUI 复现

## 功能总览

### 对接引擎
- **AutoDock Vina 1.2.7** 开箱内置，无需任何安装
- **smina** 引擎（Vinardo 打分 / off-docking 工作流）
- **gnina / Uni-Dock** GPU 引擎可选插拔（Linux + NVIDIA），运行路径自动探测
- 引擎按平台随安装包分发；外部引擎可在设置中注册接入

### 双模式工作流
- 五步向导：受体 → 配体 → 口袋 → 对接 → 分析，每步带校验与明确指引
- 内置**一键示例项目**：HIV-1 蛋白酶（1HSG）受体/配体预置，一键跑通端到端
- 批量虚拟筛选：分子库队列、暂停/取消、亲和性排行、一键 CSV 导出
- 柔性受体残基支持，报错信息说人话
- **等价命令行复现**：每次 GUI 运行都展示其真实执行的 `vina`/`smina` 命令

### 分析与质检
- 受体/口袋/构象 3D 可视化（3Dmol.js）
- **PLIP** 相互作用剖析（氢键、疏水、π 堆积、盐桥……）
- 2D 相互作用示意图、截图导出
- **PoseBusters** 构象合理性质检
- 共晶配体 RMSD 自对接校验

### 可溯源与可复现（2.0 新增）
- 每步运行记录完整**事件流（JSONL）**
- 一键导出 **`audit.zip` 审计包**：输入文件（含 SHA256）、配置、结果、环境指纹（rdkit / meeko / openbabel 版本）、引擎版本
- **`dockstudio audit` CLI** 检查复核审计包 —— 方便 agent / 脚本自动化审计
- 依赖上界锁定（如 `rdkit<2026.4`），保证跨环境可复现

### 输入数据从哪来
- 应用内**输入格式规范指引卡**：合法的受体 PDB / 配体 SDF / SMILES 长什么样
- 来源向导直达 **RCSB PDB**、**PubChem**、**ZINC**，含一步步下载指引
- SMILES 即查 + 内置**演示分子库**（10 个真实上市药物）
- 中英双语界面

## 下载安装

到 [Releases](https://github.com/2020chenjinkui-cloud/dockstudio/releases) 页取最新安装包：

| 平台 | 文件 | 说明 |
|---|---|---|
| Windows | `DockStudio_*_x64-setup.exe`（NSIS）/ `*_x64_en-US.msi` | 免管理员安装到用户目录或指定盘；需 WebView2（Win11 自带） |
| macOS（Apple Silicon） | `DockStudio_*_aarch64.dmg` | 首次启动：右键 → 打开（未签名） |
| Linux | `*.AppImage` / `*_amd64.deb` / `*.x86_64.rpm` | AppImage 需 FUSE（新发行版装 `libfuse2`） |

GPU 引擎（gnina / Uni-Dock）目前支持 **Linux + NVIDIA**；安装包全平台内置 CPU 引擎，GPU 引擎从 `PATH` / 引擎设置中自动探测。

## 快速上手

1. 安装并启动 DockStudio
2. 选**简易模式** → 点**加载示例项目（1HSG）**
3. 五步默认参数走完 —— 得到对接构象、亲和性打分、PLIP 相互作用与 PoseBusters 报告
4. 打开**导出 → 审计包**查看溯源档案，或复制**等价命令行**到终端复跑一遍

跑自己的靶点？放入受体 PDB（RCSB 下载）和配体（PubChem 的 SDF，或直接粘 SMILES）—— 每一步的指引卡都会说明接受的格式。

## Docker

```bash
docker compose up --build
# 后端 API：http://127.0.0.1:8000（健康检查 /api/health）
```

镜像为**纯后端服务模式**（API 无界面）—— 通过 REST API / SDK 消费，或作为桌面客户端的算力后端。可选 GPU worker：`docker compose --profile gpu up`，启 gnina / Uni-Dock 实例（需 NVIDIA Container Toolkit）。

## 开发

```bash
# 后端（FastAPI + 核心域 + 引擎插件，包名 dockstudio）
cd server
pip install -e ".[dev]"
uvicorn dockstudio.api.app:create_app --factory --reload

# 前端（React 18 + Vite + TypeScript + 3Dmol.js）
cd apps/web
npm install
npm run dev

# 桌面壳（Tauri 2，Rust + 后端 sidecar）
cd apps/desktop && npm install
npm run tauri dev
```

`pip install -e server` 后提供 CLI：`dockstudio`（含 `dockstudio audit`）与 `dockstudio-server`。

### 仓库结构

```
server/      Python 后端：FastAPI API + 核心域 + 引擎插件
apps/web/    React + Vite 前端（3Dmol.js 查看器）
apps/desktop/ Tauri 2 桌面壳（Rust）
resources/   随安装包分发的引擎二进制（按平台分目录，不入 git）
scripts/     打包、冒烟、发版脚本（fetch_engines.sh、pyinstaller spec……）
docker/      Dockerfile + docker-compose.yml
docs/        开发文档：API 契约、审计包规范、SDK、部署、AI agent 指南
```

## 测试与质控

- **112 个后端测试**（pytest，含真实 Vina 端到端）+ ruff 静态检查
- **97 个前端测试**
- **10 个共晶靶点金标准重对接回归**（自托管 CI 每周自动跑）
- 安装包实测：Windows NSIS 静默安装、Linux 实体机 DEB/AppImage 安装验证
- CI：GitHub Actions（win/mac 走托管 runner，Linux 构建与 GPU 相关 job 走自托管 runner）

## 已知限制

- **Boltz-2 AI 亲和性复核**：作为可选步骤支持，但需自备 `boltz` 安装与 ≥10GB 显存 GPU，真机验证仍在排队
- macOS 安装包已构建并过 CI，但尚未在实体 Mac 上验证
- gnina / Uni-Dock GPU 引擎仅 Linux + NVIDIA

## 文档

- [API 契约](docs/api-contract.md) · [审计包规范](docs/audit-schema.md) · [SDK](docs/sdk.md) · [部署](docs/deploy.md) · [AI agent 指南](docs/ai.md)

## 许可证

尚未声明许可证 —— 保留作者全部权利。如需复用代码，请先开 issue 讨论。
