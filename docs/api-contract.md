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

---

# v0.3 增量（M2 + M1 收尾）

## 口袋预测

### POST /api/pockets
Body `{receptor_id, top_n?: 5}` → `PocketResponse {receptor_id, method:"builtin_grid", pockets:[{id, center, size, score, n_points}]}`。
内置几何算法（网格埋藏腔检测），几秒级；结果可缓存（receptor_id+参数）。失败 400。

## 柔性残基

### POST /api/prepare/receptor/{rid}/flex
Body `{flex_residues: ["A:42", ...]}`（`链:残基号`，<=10 个）→ 新的 `PreparedReceptor` 变体
（parent_id=rid，flex_residues、flex_pdbqt_path、receptor_json_path 非空）。残基不存在 → 400 并列出可用链。
柔性变体参与对接时，引擎须带 `--flex`；不支持 flex 的引擎用柔性变体会 400。

## 配体库

### POST /api/prepare/library
multipart `file`（.sdf 多分子 / .smi / .csv 含 smiles 列），可选 `name` → `PreparedLibrary {id, name, path, ligand_ids[], failed:[{index,smiles,reason}]}`。
流式解析，单条失败不中断；上限 50000 条。`.smi` 每行 `SMILES [名称]`。

## 批量筛选

### POST /api/projects/{pid}/batches
Body = `BatchCreate {name?, engine, receptor_id, ligand_ids[], box, params}` → `{batch_id}`。
为每个配体建 Job（batch_id 标记）并入同一串行队列。
### GET /api/batches/{bid} → `BatchSummary`（由 jobs 聚合：done/failed/running/queued/canceled、is_terminal）
### POST /api/batches/{bid}/cancel → 取消该批全部非终态 job（queued 直接终态，running 尽力杀）
### GET /api/batches/{bid}/results?sort=affinity&order=asc&limit=50&offset=0
→ `{total, rows:[BatchRow{ligand_id, ligand_name, smiles, job_id, status, top_affinity, error}]}`
### GET /api/batches/{bid}/results.csv → text/csv（列：rank, ligand_name, smiles, top_affinity, status, job_id, error）
### WS /ws/batches/{bid} → `{"type":"batch_progress","data":BatchSummary}`，终态发 `{"type":"finished","data":BatchSummary}`

## 其他

- GET /api/jobs/{jid} 响应 Job 含 batch_id 字段（可空）
- GET /api/prepare/receptors 返回含 flex 变体（前端应区分显示 parent_id 非空者）
- CLI：`py -3.12 -m dockstudio.cli serve|engines|demo`（pyproject 增加 script `dockstudio`）

---

# v0.4 增量（M3）

## 姿态质检与 2D 相互作用图

### GET /api/jobs/{jid}/analysis?pose_index=N
响应 PoseAnalysis 新增 `pose_check: PoseCheck|null.dtype|null`
（`{pose_index, ok, source:"posebusters"|"builtin", issues:[{rule, severity, message}]}`）。

### GET /api/jobs/{jid}/diagram?pose_index=N
→ `image/svg+xml`：2D 相互作用示意图（配体 2D + 周围残基标签，按相互作用类型上色）。

## 报告 v2

### GET /api/jobs/{jid}/report?sections=snapshot,scores,interactions,diagram,viewer3d,disclaimer&title=...
- sections：逗号分隔子集（默认全部）；title：自定义中文标题（<=80 字符，缺省用项目/任务名）。

## 结果对比

### GET /api/compare?job_ids=id1,id2,id3
→ `ComparisonResult {rows:[ComparisonRow{job_id, ligand_id, ligand_name, smiles, engine, scoring, exhaustiveness, top_affinity, n_poses, status}]}`。
跨项目可混选；不存在的 id 跳过并在 rows 中不出现（总数 mismatch 即知）。

## 统计

### GET /api/stats
→ `{projects:int, jobs:int, batches:int, completed_jobs:int, avg_dock_seconds:float|null, by_engine:{vina:int,...}}`。

## 引擎（/api/engines 体现，无需新端点）
- 新增 gnina 适配器：manifest platforms=["linux"]，requires_gpu=true（CPU 可跑但极慢，unavailable_reason 如实说明，WSL2 用户按指引自查）；scoring_functions=["vina","vinardo","ad4_scoring","dkoes_fast"]，supports_flex_residues。
- 新增 unidock 适配器：platforms=["linux"]，requires_gpu=true，supports_batch_ligands=true。
- 二进制解析顺序统一：env(var) → PATH → resources/bin/<platform>/。
