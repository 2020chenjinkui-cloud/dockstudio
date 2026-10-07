import { useEffect, useState } from 'react';
import { getProject } from '../api';
import { statusLabel } from '../labels';
import type { Project, ProjectDetail, ProjectJobSummary } from '../types';
import ConfirmDialog from './ConfirmDialog';

interface ProjectBarProps {
  projects: Project[];
  current: Project | null;
  backendOk: boolean;
  onSelect: (id: string) => void;
  onCreate: (name: string) => Promise<void>;
  onDelete: () => Promise<void>;
  onOpenJob: (jobId: string) => void;
  onOpenBatch: (batchId: string) => void;
}

type DrawerItem =
  | { kind: 'job'; job: ProjectJobSummary }
  | { kind: 'batch'; batchId: string; jobs: ProjectJobSummary[] };

function groupJobs(jobs: ProjectJobSummary[]): DrawerItem[] {
  const items: DrawerItem[] = [];
  const batches = new Map<string, { index: number; jobs: ProjectJobSummary[] }>();
  for (const job of jobs) {
    if (job.batch_id) {
      const entry = batches.get(job.batch_id);
      if (entry) {
        entry.jobs.push(job);
      } else {
        batches.set(job.batch_id, { index: items.length, jobs: [job] });
        items.push({ kind: 'batch', batchId: job.batch_id, jobs: [job] });
      }
    } else {
      items.push({ kind: 'job', job });
    }
  }
  return items;
}

function batchTopAffinity(jobs: ProjectJobSummary[]): number | null {
  const values = jobs.map((j) => j.top_affinity).filter((v): v is number => v !== null);
  if (values.length === 0) return null;
  return Math.min(...values);
}

export default function ProjectBar({
  projects,
  current,
  backendOk,
  onSelect,
  onCreate,
  onDelete,
  onOpenJob,
  onOpenBatch,
}: ProjectBarProps) {
  const [creating, setCreating] = useState(false);
  const [newName, setNewName] = useState('');
  const [createError, setCreateError] = useState<string | null>(null);
  const [drawerOpen, setDrawerOpen] = useState(false);
  const [detail, setDetail] = useState<ProjectDetail | null>(null);
  const [confirmDelete, setConfirmDelete] = useState(false);

  useEffect(() => {
    if (!drawerOpen || !current) return;
    let cancelled = false;
    getProject(current.id)
      .then((d) => {
        if (!cancelled) setDetail(d);
      })
      .catch(() => {
        if (!cancelled) setDetail(null);
      });
    return () => {
      cancelled = true;
    };
  }, [drawerOpen, current]);

  const submitCreate = async () => {
    const name = newName.trim();
    if (!name) return;
    setCreateError(null);
    try {
      await onCreate(name);
      setNewName('');
      setCreating(false);
    } catch (e) {
      setCreateError(e instanceof Error ? e.message : String(e));
    }
  };

  const items = detail ? groupJobs(detail.jobs) : [];

  return (
    <header className="project-bar">
      <span className="brand">DockStudio</span>
      {current ? (
        <div className="project-controls">
          <select value={current.id} onChange={(e) => onSelect(e.target.value)}>
            {projects.map((p) => (
              <option key={p.id} value={p.id}>
                {p.name}
              </option>
            ))}
          </select>
          {creating ? (
            <span className="create-form">
              <input
                type="text"
                placeholder="新项目名称"
                value={newName}
                autoFocus
                onChange={(e) => setNewName(e.target.value)}
                onKeyDown={(e) => {
                  if (e.key === 'Enter') void submitCreate();
                  if (e.key === 'Escape') setCreating(false);
                }}
              />
              <button onClick={() => void submitCreate()}>创建</button>
              <button onClick={() => setCreating(false)}>取消</button>
            </span>
          ) : (
            <button onClick={() => setCreating(true)}>新建项目</button>
          )}
          <button onClick={() => setDrawerOpen(true)}>任务历史</button>
          <button className="danger-outline" onClick={() => setConfirmDelete(true)}>
            删除项目
          </button>
        </div>
      ) : (
        <span className="muted">未选择项目</span>
      )}
      <span className={backendOk ? 'backend-status ok' : 'backend-status bad'}>
        {backendOk ? '后端已连接' : '后端未启动'}
      </span>
      {createError && <span className="error-text">{createError}</span>}

      {drawerOpen && (
        <div className="drawer-overlay" onClick={() => setDrawerOpen(false)}>
          <aside className="drawer" onClick={(e) => e.stopPropagation()}>
            <div className="drawer-head">
              <h3>任务历史</h3>
              <button className="link-btn" onClick={() => setDrawerOpen(false)}>
                关闭
              </button>
            </div>
            {!detail && <p className="muted">加载中…</p>}
            {detail && items.length === 0 && <p className="muted">本项目还没有任务。</p>}
            {detail && (
              <ul className="job-list">
                {[...items].reverse().map((item) =>
                  item.kind === 'job' ? (
                    <li key={item.job.id}>
                      <button
                        className="job-item"
                        onClick={() => {
                          onOpenJob(item.job.id);
                          setDrawerOpen(false);
                        }}
                      >
                        <span className={`badge status-${item.job.status}`}>
                          {statusLabel(item.job.status)}
                        </span>
                        <span className="job-meta">
                          <span className="mono">{item.job.id.slice(0, 8)}</span>
                          <span className="muted">{item.job.engine}</span>
                        </span>
                        <span className="job-affinity">
                          {item.job.top_affinity !== null
                            ? `${item.job.top_affinity.toFixed(1)} kcal/mol`
                            : '—'}
                        </span>
                      </button>
                    </li>
                  ) : (
                    <li key={item.batchId}>
                      <button
                        className="job-item"
                        onClick={() => {
                          onOpenBatch(item.batchId);
                          setDrawerOpen(false);
                        }}
                      >
                        <span className="badge badge-flex">批量</span>
                        <span className="job-meta">
                          <span className="mono">{item.batchId.slice(0, 8)}</span>
                          <span className="muted">{item.jobs.length} 配体</span>
                        </span>
                        <span className="job-affinity">
                          {batchTopAffinity(item.jobs) !== null
                            ? `top ${batchTopAffinity(item.jobs)!.toFixed(1)}`
                            : '—'}
                        </span>
                      </button>
                    </li>
                  ),
                )}
              </ul>
            )}
          </aside>
        </div>
      )}

      {confirmDelete && current && (
        <ConfirmDialog
          title="删除项目"
          message={`确定删除项目「${current.name}」吗？项目将被标记删除，任务历史不再可见。`}
          confirmText="删除"
          onCancel={() => setConfirmDelete(false)}
          onConfirm={() => {
            setConfirmDelete(false);
            void onDelete();
          }}
        />
      )}
    </header>
  );
}
