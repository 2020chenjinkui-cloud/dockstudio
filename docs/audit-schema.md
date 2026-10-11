# DockStudio 事件 Schema（events.jsonl）与 audit.zip 结构

每个项目根目录一条 `events.jsonl`（数据目录下 `projects/<pid>/events.jsonl`；受体/配体
准备产生的事件记在伪项目 `_prepared` 下）。每行一个 JSON 对象：

```json
{"ts":"2026-10-11T12:34:56.789+00:00","actor":"user","action":"job.created","subject":{"kind":"job","id":"a1b2c3d4"},"payload":{"engine":"vina","box":{"center":[13.0,22.4,5.5],"size":[22.5,22.5,22.5]},"params":{"exhaustiveness":4},"env":{"python":"3.12.0","rdkit":"2026.03.6","meeko":"0.8.0","openbabel":"3.1.1","platform":"Windows"},"input_sha256":{"receptor":"...","ligand":"..."}}}
```

| 字段 | 说明 |
|---|---|
| `ts` | ISO8601 UTC（毫秒），附加一次性 |
| `actor` | `user`（API 动作）/ `system`（队列/worker 生命周期） |
| `action` | 见下表动作全集 |
| `subject` | `{kind, id}`，产生对象 |
| `payload` | 动作相关全量数据，用于 agent 审计 |

## 动作全集

| action | subject | payload 关键字段 |
|---|---|---|
| `project.create` | project | `name` |
| `project.delete` | project | `name`（记录在 .trash- 副本） |
| `receptor.prepared` | receptor | `filename, original_sha256, n_chains, chains, n_cocrystals, cocrystals[], warnings[]` |
| `receptor.flex_variant` | receptor | `parent_id, flex_residues[]` |
| `ligand.prepared` | ligand | `source, smiles, n_rotatable, name` |
| `library.prepared` | library | `name, n_ligands, n_failed, filename` |
| `job.created` | job | `engine, receptor_id, ligand_id, box, params, batch_id, env(EnvFingerprint), input_sha256{receptor,ligand}` |
| `job.started` | job | `engine`（actor=system） |
| `job.completed` | job | `top_affinity, n_poses`（actor=system） |
| `job.failed` | job | `error`（用户可读中文，actor=system） |
| `job.canceled` | job | 空（actor=system） |
| `batch.created` | batch | `name, total, receptor_id, engine, ligand_ids[]` |
| `batch.canceled` | batch | 空 |
| `batch.completed` | batch | （终态聚合，未单独发，使用 GET /api/batches 查询） |
| `rescore.created` | rescore | `job_id, pose_indices[], msa_server_used` |
| `rescore.completed` | rescore | `n_results, results[]` |
| `rescore.canceled` | rescore | 空 |

文件层保证：同路径 append 走加锁 + O_APPEND，恶意并发行可能被截断（崩溃残留），
消费端必须容忍非 JSON 行（`read()` 自动跳过）。

## audit.zip 结构

`GET /api/jobs/{jid}/audit.zip`，文件名 `audit-<jid>.zip`：

```
audit.zip
├── snapshot.json     # {"job": <job.json 全量，含 env>,"inputs":{"receptor": "<sha256>","ligand":"<sha256>"}}
├── events.jsonl      # 该 job（及其项目）相关事件，按时间顺序
├── engine.log        # 引擎 stdout+stderr 原文
├── receptor.pdbqt    # 进 run 目录中的受体输入（存在才装）
├── ligand.pdbqt      # 配体输入（存在才装）
├── out.pdbqt         # vina 输出全部姿态（存在才装）
└── report.md         # 事理链 Markdown（见下节）
```

## report.md 内容结构（发给 agent 复核）

1. **起始终起**：项目/任务/创建-开始-结束时间/状态/错误。
2. **输入**：receptor_id/ligand_id、各输入 SHA256（如其前缀）。
3. **环境指纹**：python/rdkit/meeko/openbabel/platform（复现对照）。
4. **盒子与参数**：center/size、engine/scoring/exhaustiveness/num_modes/energy_range/cpu/seed。
5. **姿态打分**：全部 pose 的 affinity 与 RMSD 区间。
6. **事件序列**：每条事件的 ts/action/actor/payload 摘要。
7. **免责声明**：打分仅作相对排序，AI 复核为独立预测，不代表实测亲和力。

## 审计建议

- 复现对照：先比 `env` 与 `input_sha256`；不一致=输入/环境已变。
- 争议姿态：用 `GET /api/jobs/{jid}/crystal_rmsd?pose_index=N` 与共晶构象做正交校验
  （重对接 rmsd<2Å 视为高质量，2-4Å 可接受，>4Å 警惕）。
- 侧信道：当 `msa_server_used=true`，上游蛋白序列已发送至公共 ColabFold MSA 服务器。
