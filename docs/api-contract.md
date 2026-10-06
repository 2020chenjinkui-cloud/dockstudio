# API 契约 v0.2（M1）

前后端以此文件为准。错误统一 `{"detail": "<用户可读信息>"}`。
M0 端点（health / engines / prepare / dock / file）保持不变，见文末；M1 新增项目、异步任务、分析、报告。

## 项目 Projects

### POST /api/projects
Body `{name: string}` → `Project {id, name, path, created_at}`。
项目目录 = data_dir/projects/<id>/，内含 jobs/<job_id>/。

### GET /api/projects
→ `[Project, ...]`（按创建时间倒序）

### GET /api/projects/{pid}
→ Project + `jobs: [JobSummary]`。JobSummary = `{id, status, ligand_id, engine, top_affinity: number|null, created_at, finished_at}`。

### DELETE /api/projects/{pid}
软删（目录重命名加 `.trash-` 前缀）。204。

## 准备物列表（向导选择器用）

### GET /api/prepare/receptors → `[PreparedReceptor]`
### GET /api/prepare/ligands → `[PreparedLigand]`

## 异步任务 Jobs

### POST /api/projects/{pid}/jobs
Body = `DockJobRequest {engine, receptor_id, ligand_id, box, params}` → `{job_id}`。
立即入队。队列串行执行（M1），重启后未完成任务标记 interrupted。

### GET /api/jobs/{jid}
→ `Job {id, project_id, status, progress:{stage,percent,message}, request, result:DockResult|null, error, created_at, started_at, finished_at}`

### POST /api/jobs/{jid}/cancel
queued → 直接取消；running → 调引擎 cancel_current() 尽力杀进程。→ Job。

### WS /ws/jobs/{jid}
推送 JSON 事件：`{"type":"status","data":<Job>}`、`{"type":"progress","data":<JobProgress>}`、
`{"type":"finished","data":<Job>}`。客户端断线后应回到 GET /api/jobs/{jid} 轮询。

## 结果分析

### GET /api/jobs/{jid}/analysis?pose_index=1
→ `PoseAnalysis {pose_index, source:"plip"|"builtin", interactions:[{type, residue, distance, ligand_atom_ids, detail}], warnings}`。
优先 PLIP，不可用则内置几何检测（氢键/疏水/π堆积，距离+角度规则），结果缓存在 job 目录。

### GET /api/jobs/{jid}/report
→ `text/html` 自包含报告（内嵌参数快照、打分表、相互作用、3D 场景（3Dmol CDN + 内联结构数据）、免责声明与引用）。Content-Disposition: attachment。

## M0 保留端点（不变，向后兼容调试页）

- GET /api/health → `{status, version}`
- GET /api/engines → `[manifest + available + unavailable_reason]`
- POST /api/prepare/receptor（multipart file）→ PreparedReceptor
- POST /api/prepare/ligand（multipart file 或 smiles/name）→ PreparedLigand
- POST /api/dock（同步，DockJobRequest 同构）→ DockResult（M0 调试页继续用；不入队）
- GET /api/file?path=（data_dir 沙箱）
