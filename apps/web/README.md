# DockStudio Web

React 18 + Vite + TypeScript 前端。契约见 `docs/api-contract.md`（当前 v0.5）。

## 路由

- `/` — 主界面：五步向导（受体 → 配体 → 对接位点 → 运行 → 结果）
  - 顶部项目栏：切换/新建/删除项目、任务历史抽屉（batch 按 `batch_id` 分组为"批量·N 配体·top X"条目，可回到批量进度/排行榜）、/api/stats 统计小卡、中英文切换（localStorage 记忆）
  - 模式切换：单一对接 / 批量筛选（配体步随之变为库上传或多选）
  - 高级模式开关：完整参数（energy_range、seed）+ 柔性残基卡（`链:残基号` 客户端校验，提交前自动创建 flex 变体）
  - 对接位点步：口袋预测（/api/pockets，按 score 排序 + 横条可视化，点选定盒）
  - 任务进度：WS 订阅，断线自动降级 2s 轮询（job 与 batch 同理）
  - 结果页：姿态表联动 3D、姿态质检徽章（pose_check，最高 severity 配色，点击弹出 issues，"只看合格"过滤）、相互作用残基高亮与 2D 示意图 tab（/api/jobs/{jid}/diagram）、等效命令行卡（可复制）、报告 v2 下载对话框（章节勾选 + 自定义标题）
  - AI 复核（Boltz-2）：按 /api/capabilities 置灰（tooltip 安装指引 + docs/ai.md 提示）；对话框多选姿态（默认 top 3）+ MSA 服务器开关（默认开，红字提示序列外发）；3s 轮询进度 + 可取消；姿态表合并最新复核结果（logIC50 三位小数 / binder% 一位，表头与单元格 hover 注明 AI 独立预测）；复核历史列表可取消进行中任务；MSA 使用记录常驻红字提示
  - 结果对比：历史抽屉勾选 ≥2 个已完成任务 → 对比模态（最优 affinity 加粗高亮，前端导出 CSV，跨项目可混选）
  - 批量进度页：分段进度条（done/running/failed/canceled/queued）、排行榜（分页 + affinity 排序 + 导 CSV + 行跳 job 结果页）
- `/debug` — M0 调试台（prepare/dock/file 原始链路，保留供排查）

## i18n（zh/en）

`src/i18n.ts` 自研轻量方案：字典 + `translate()` + `I18nContext`/`useI18n()`，语言存
localStorage（键 `dockstudio.lang`），默认中文。覆盖：向导步骤名、模式/高级切换、导航按钮、
状态徽章、阶段映射、项目栏、帮助 tooltip、姿态质检文案、报告对话框、对比表头与按钮。
未覆盖（保持中文）：免责声明长文案、onboarding 正文、/debug 调试台、后端下发的错误 detail。

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

- `src/types.ts` — 与 API 契约（v0.4）一一对应的类型定义
- `src/api.ts` — typed client（health/engines/prepare/flex/dock/file/projects/jobs/batches/pockets/library/analysis/diagram/report/compare/stats）+ `subscribeJob` / `subscribeBatch`（WS 优先、断线轮询）
- `src/i18n.ts` — zh/en 字典、`translate()`、React context
- `src/labels.ts` — 状态/阶段/相互作用中文映射、残基标签解析、质检徽章映射
- `src/flex.ts` — 柔性残基输入解析与校验（格式、链存在性、去重、上限 10）
- `src/cli.ts` — 由 DockRequest 渲染等效 vina 命令行
- `src/compare.ts` — 对比选择逻辑与 CSV 生成
- `src/rescore.ts` — AI 复核结果合并（最新一次完成任务为权威）、默认姿态选择、表单校验
- `src/helpText.ts` — 帮助与提示文案（zh 单一事实来源，i18n zh 包引用）
- `src/router.ts` — 极简 history 路由
- `src/MolViewer.tsx` — 3Dmol.js 查看器（cartoon + 口袋 stick + 半透明盒、姿态/残基高亮）
- `src/App.tsx` — I18nProvider、路由分发与主界面壳（项目选择、空状态、onboarding、对比模态）
- `src/DebugPage.tsx` — M0 调试台
- `src/components/` — ProjectBar、ConfirmDialog、OnboardingOverlay、Tooltip、CliCard、StatsChip、ReportDialog、CompareDialog、AiRescoreSection、AiRescoreDialog 及 `wizard/`（Wizard、五个步骤、StepBatch、BatchPanel）
- `src/*.test.ts` — vitest 单测

## TODO

- TODO(M1 遗留): 3D 查看器由 3Dmol.js 迁移到 Mol*；`/api/file` 调试接口替换为项目内文件接口。
- TODO: `GET /api/prepare/receptors` 若一次返回大量 flex 变体，库列表需加过滤（变体/原始分组）。
- TODO: 排行榜排序目前固定 affinity asc，等后端确认更多 sort 字段后再开放表头排序。
- TODO: onboarding 与 /debug 页尚未接 i18n；免责声明长文案待后续统一走文案资源。
- TODO(M4 联调点): `binder_probability` 前端按 0–1 处理（×100 显示百分比）；若后端直接给 0–100 需去掉换算。
- TODO(M4 联调点): 复核合并视图取"最新一次有结果的任务"（按 finished_at，缺失回退 created_at）；若后端采用"最后一次成功"语义则一致，若需跨任务合并多姿态覆盖再调 `mergeRescoreResults`。
- TODO: 复核进度轮询间隔 3s（裸 HTTP，无 WS）；v0.5 契约未给复核 WS，后续若加可复用 subscribeJob 模式。
