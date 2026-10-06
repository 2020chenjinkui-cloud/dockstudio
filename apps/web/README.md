# DockStudio Web

React 18 + Vite + TypeScript 前端。契约见 `docs/api-contract.md`（当前 v0.2）。

## 路由

- `/` — 简易模式：五步向导（受体 → 配体 → 对接位点 → 运行 → 结果）+ 顶部项目栏
  （切换/新建/删除项目、任务历史抽屉）+ 任务进度（WS 订阅，断线自动降级 2s 轮询）+
  结果页（姿态表、3D 相互作用残基高亮、HTML 报告下载）+ 首次使用三步引导
- `/debug` — M0 调试台（prepare/dock/file 原始链路，保留供排查）

## 启动开发环境

```bash
# 1. 先启动后端（约定 http://localhost:8000）
cd server
pip install -e ".[dev]"
uvicorn dockstudio.api.app:create_app --factory --reload

# 2. 启动前端（另开终端，/api 与 /ws 由 Vite 代理到 8000 端口）
cd apps/web
npm install
npm run dev   # http://localhost:5173
```

后端未启动时页面仍可打开：项目区显示"后端未启动"，各请求错误以横幅展示 detail 或
"后端未启动或无法连接"。

## 常用命令

- `npm run dev` — 开发服务器
- `npm run build` — `tsc --noEmit` + `vite build`
- `npm run lint` — `tsc --noEmit`（未引入 eslint）
- `npm test` — vitest（api client 与状态/标签映射单测）
- `npm run preview` — 预览构建产物

## 代码结构

- `src/types.ts` — 与 API 契约（v0.2）一一对应的类型定义
- `src/api.ts` — typed client（health/engines/prepare/dock/file/projects/jobs/analysis/report）+ `subscribeJob`（WS 优先、断线轮询）
- `src/labels.ts` — 任务状态/阶段/相互作用类型的中文映射、残基标签解析（`ASP25:A`）
- `src/helpText.ts` — 参数帮助与提示文案（集中维护）
- `src/router.ts` — 极简 history 路由
- `src/MolViewer.tsx` — 3Dmol.js 查看器（受体 cartoon + 口袋 stick + 半透明盒、姿态高亮、残基高亮）
- `src/App.tsx` — 路由分发与向导页壳（项目选择、空状态、onboarding）
- `src/DebugPage.tsx` — M0 调试台
- `src/components/` — ProjectBar、ConfirmDialog、OnboardingOverlay、Tooltip 及 `wizard/` 五个步骤组件
- `src/*.test.ts` — vitest 单测

## TODO

- TODO(M1 遗留): 3D 查看器由 3Dmol.js 迁移到 Mol*；`/api/file` 调试接口替换为项目内文件接口。
- TODO(M2): "高级模式"入口（目前为 disabled 占位）。
- TODO(M3): "按口袋建议"定盒（自动口袋预测）。
