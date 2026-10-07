# DockStudio Web

React 18 + Vite + TypeScript 前端。契约见 `docs/api-contract.md`（当前 v0.3）。

## 路由

- `/` — 主界面：五步向导（受体 → 配体 → 对接位点 → 运行 → 结果）
  - 顶部项目栏：切换/新建/删除项目、任务历史抽屉（batch 按 `batch_id` 分组为"批量·N 配体·top X"条目，可回到批量进度/排行榜）
  - 模式切换：单一对接 / 批量筛选（配体步随之变为库上传或多选）
  - 高级模式开关：完整参数（energy_range、seed）+ 柔性残基卡（`链:残基号` 客户端校验，提交前自动创建 flex 变体）
  - 对接位点步：口袋预测（/api/pockets，按 score 排序 + 横条可视化，点选定盒）
  - 任务进度：WS 订阅，断线自动降级 2s 轮询（job 与 batch 同理）
  - 结果页：姿态表联动 3D、相互作用残基高亮、等效命令行卡（可复制）、HTML 报告下载
  - 批量进度页：分段进度条（done/running/failed/canceled/queued）、排行榜（分页 + affinity 排序 + 导 CSV + 行跳 job 结果页）
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
- `npm test` — vitest（api client、标签映射、flex 校验、CLI 渲染）
- `npm run preview` — 预览构建产物

## 代码结构

- `src/types.ts` — 与 API 契约（v0.3）一一对应的类型定义
- `src/api.ts` — typed client（health/engines/prepare/flex/dock/file/projects/jobs/batches/pockets/library/analysis/report）+ `subscribeJob` / `subscribeBatch`（WS 优先、断线轮询）
- `src/labels.ts` — 状态/阶段/相互作用中文映射、残基标签解析
- `src/flex.ts` — 柔性残基输入解析与校验（格式、链存在性、去重、上限 10）
- `src/cli.ts` — 由 DockRequest 渲染等效 vina 命令行
- `src/helpText.ts` — 参数帮助与提示文案（集中维护）
- `src/router.ts` — 极简 history 路由
- `src/MolViewer.tsx` — 3Dmol.js 查看器（cartoon + 口袋 stick + 半透明盒、姿态/残基高亮）
- `src/App.tsx` — 路由分发与主界面壳（项目选择、空状态、onboarding）
- `src/DebugPage.tsx` — M0 调试台
- `src/components/` — ProjectBar、ConfirmDialog、OnboardingOverlay、Tooltip、CliCard 及 `wizard/`（Wizard、五个步骤、StepBatch、BatchPanel）
- `src/*.test.ts` — vitest 单测

## TODO

- TODO(M1 遗留): 3D 查看器由 3Dmol.js 迁移到 Mol*；`/api/file` 调试接口替换为项目内文件接口。
- TODO: `GET /api/prepare/receptors` 若一次返回大量 flex 变体，库列表需加过滤（变体/原始分组）。
- TODO: 排行榜排序目前固定 affinity asc，等后端确认更多 sort 字段后再开放表头排序。
