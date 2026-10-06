import { useEffect, useState } from 'react';
import { getProject } from '../api';
import { statusLabel } from '../labels';
import type { Project, ProjectDetail } from '../types';
import ConfirmDialog from './ConfirmDialog';

interface ProjectBarProps {
  projects: Project[];
  current: Project | null;
  backendOk: boolean;
  onSelect: (id: string) => void;
  onCreate: (name: string) => Promise<void>;
  onDelete: () => Promise<void>;
  onOpenJob: (jobId: string) => void;
}

export default function ProjectBar({
  projects,
  current,
  backendOk,
  onSelect,
  onCreate,
  onDelete,
  onOpenJob,
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
            {detail && detail.jobs.length === 0 && <p className="muted">本项目还没有任务。</p>}
            {detail && (
              <ul className="job-list">
                {[...detail.jobs].reverse().map((j) => (
                  <li key={j.id}>
                    <button
                      className="job-item"
                      onClick={() => {
                        onOpenJob(j.id);
                        setDrawerOpen(false);
                      }}
                    >
                      <span className={`badge status-${j.status}`}>{statusLabel(j.status)}</span>
                      <span className="job-meta">
                        <span className="mono">{j.id.slice(0, 8)}</span>
                        <span className="muted">{j.engine}</span>
                      </span>
                      <span className="job-affinity">
                        {j.top_affinity !== null ? `${j.top_affinity.toFixed(1)} kcal/mol` : '—'}
                      </span>
                    </button>
                  </li>
                ))}
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
