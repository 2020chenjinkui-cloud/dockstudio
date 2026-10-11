import { useEffect, useState } from 'react';
import { getProject, projectEventsUrl } from '../api';
import { toggleCompareId } from '../compare';
import { useI18n } from '../i18n';
import type { Project, ProjectDetail, ProjectJobSummary } from '../types';
import ConfirmDialog from './ConfirmDialog';
import StatsChip from './StatsChip';

interface ProjectBarProps {
  projects: Project[];
  current: Project | null;
  backendOk: boolean;
  compareIds: string[];
  onCompareIdsChange: (ids: string[]) => void;
  onOpenCompare: () => void;
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
  const batches = new Map<string, { jobs: ProjectJobSummary[] }>();
  for (const job of jobs) {
    if (job.batch_id) {
      const entry = batches.get(job.batch_id);
      if (entry) {
        entry.jobs.push(job);
      } else {
        batches.set(job.batch_id, { jobs: [job] });
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
  compareIds,
  onCompareIdsChange,
  onOpenCompare,
  onSelect,
  onCreate,
  onDelete,
  onOpenJob,
  onOpenBatch,
}: ProjectBarProps) {
  const { lang, setLang, t } = useI18n();
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
                placeholder={t('bar.namePlaceholder')}
                value={newName}
                autoFocus
                onChange={(e) => setNewName(e.target.value)}
                onKeyDown={(e) => {
                  if (e.key === 'Enter') void submitCreate();
                  if (e.key === 'Escape') setCreating(false);
                }}
              />
              <button onClick={() => void submitCreate()}>{t('bar.create')}</button>
              <button onClick={() => setCreating(false)}>{t('bar.cancel')}</button>
            </span>
          ) : (
            <button onClick={() => setCreating(true)}>{t('bar.new')}</button>
          )}
          <button onClick={() => setDrawerOpen(true)}>{t('bar.history')}</button>
          <button className="danger-outline" onClick={() => setConfirmDelete(true)}>
            {t('bar.delete')}
          </button>
        </div>
      ) : (
        <span className="muted">{t('bar.none')}</span>
      )}
      <span className="bar-right">
        <StatsChip />
        <span className={backendOk ? 'backend-status ok' : 'backend-status bad'}>
          {backendOk ? t('bar.backendOk') : t('bar.backendDown')}
        </span>
        <span className="lang-switch">
          <button
            className={lang === 'zh' ? 'chip active' : 'chip'}
            onClick={() => setLang('zh')}
          >
            中
          </button>
          <button
            className={lang === 'en' ? 'chip active' : 'chip'}
            onClick={() => setLang('en')}
          >
            EN
          </button>
        </span>
      </span>
      {createError && <span className="error-text">{createError}</span>}

      {drawerOpen && (
        <div className="drawer-overlay" onClick={() => setDrawerOpen(false)}>
          <aside className="drawer" onClick={(e) => e.stopPropagation()}>
            <div className="drawer-head">
              <h3>{t('bar.history')}</h3>
              <span className="drawer-head-actions">
                {current && (
                  <a
                    className="link-btn"
                    href={projectEventsUrl(current.id)}
                    target="_blank"
                    rel="noreferrer"
                  >
                    {t('audit.events')}
                  </a>
                )}
                <button className="link-btn" onClick={() => setDrawerOpen(false)}>
                  {t('bar.close')}
                </button>
              </span>
            </div>
            <div className="drawer-body">
              {!detail && <p className="muted">{t('common.loading')}</p>}
              {detail && items.length === 0 && <p className="muted">{t('bar.jobsEmpty')}</p>}
              {detail && (
                <ul className="job-list">
                  {[...items].reverse().map((item) =>
                    item.kind === 'job' ? (
                      <li key={item.job.id} className="job-row">
                        {item.job.status === 'completed' && (
                          <input
                            type="checkbox"
                            className="cmp-check"
                            title={t('cmp.title')}
                            checked={compareIds.includes(item.job.id)}
                            onChange={() =>
                              onCompareIdsChange(toggleCompareId(compareIds, item.job.id))
                            }
                          />
                        )}
                        <button
                          className="job-item"
                          onClick={() => {
                            onOpenJob(item.job.id);
                            setDrawerOpen(false);
                          }}
                        >
                          <span className={`badge status-${item.job.status}`}>
                            {t(`status.${item.job.status}`, item.job.status)}
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
                      <li key={item.batchId} className="job-row">
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
            </div>
            {compareIds.length >= 2 && (
              <div className="drawer-foot">
                <button
                  className="primary cmp-button"
                  onClick={() => {
                    onOpenCompare();
                    setDrawerOpen(false);
                  }}
                >
                  {t('cmp.button')}（{compareIds.length}）
                </button>
              </div>
            )}
          </aside>
        </div>
      )}

      {confirmDelete && current && (
        <ConfirmDialog
          title={t('bar.delete')}
          message={`确定删除项目「${current.name}」吗？项目将被标记删除，任务历史不再可见。`}
          confirmText={t('bar.delete')}
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
