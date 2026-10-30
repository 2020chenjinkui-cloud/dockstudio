import { useEffect, useRef, useState } from 'react';
import {
  cancelRescoreTask,
  getCapabilities,
  getRescoreTask,
  listRescoreTasks,
} from '../api';
import { useI18n } from '../i18n';
import { isTerminalStatus } from '../labels';
import { latestTaskWithResults } from '../rescore';
import type { AiAvailability, AiRescoreTask } from '../types';
import AiRescoreDialog from './AiRescoreDialog';

interface AiRescoreSectionProps {
  jobId: string;
  poseCount: number;
  onTasks: (tasks: AiRescoreTask[]) => void;
}

const POLL_MS = 3000;

export default function AiRescoreSection({ jobId, poseCount, onTasks }: AiRescoreSectionProps) {
  const { t } = useI18n();
  const [availability, setAvailability] = useState<AiAvailability | null>(null);
  const [tasks, setTasks] = useState<AiRescoreTask[]>([]);
  const [dialogOpen, setDialogOpen] = useState(false);
  const [runningTask, setRunningTask] = useState<AiRescoreTask | null>(null);
  const [pollError, setPollError] = useState<string | null>(null);
  const tasksRef = useRef(tasks);
  tasksRef.current = tasks;

  const pushTasks = (next: AiRescoreTask[]) => {
    setTasks(next);
    onTasks(next);
  };

  const refreshTasks = async () => {
    try {
      const list = await listRescoreTasks(jobId);
      pushTasks(list);
      const active = list.find((x) => !isTerminalStatus(x.status));
      if (active) {
        setRunningTask((cur) => (cur && !isTerminalStatus(cur.status) ? cur : active));
      }
    } catch {
      // 列表失败不阻塞页面
    }
  };

  useEffect(() => {
    getCapabilities()
      .then((c) => setAvailability(c.ai_rescore))
      .catch(() => setAvailability(null));
    void refreshTasks();
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [jobId]);

  useEffect(() => {
    if (!runningTask || isTerminalStatus(runningTask.status)) return;
    let cancelled = false;
    let timer: number | null = null;
    const tick = async () => {
      try {
        const task = await getRescoreTask(runningTask.id);
        if (cancelled) return;
        setRunningTask(task);
        if (isTerminalStatus(task.status)) {
          await refreshTasks();
          return;
        }
      } catch (e) {
        if (!cancelled) setPollError(e instanceof Error ? e.message : String(e));
      }
      if (!cancelled) timer = window.setTimeout(tick, POLL_MS);
    };
    timer = window.setTimeout(tick, POLL_MS);
    return () => {
      cancelled = true;
      if (timer !== null) clearTimeout(timer);
    };
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [runningTask?.id, runningTask?.status, jobId]);

  const cancelTask = async (target: AiRescoreTask) => {
    try {
      const task = await cancelRescoreTask(target.id);
      setRunningTask((cur) => (cur?.id === task.id ? task : cur));
      await refreshTasks();
    } catch (e) {
      setPollError(e instanceof Error ? e.message : String(e));
    }
  };

  const latestWithResults = latestTaskWithResults(tasks);
  const aiAvailable = availability?.available === true;

  return (
    <div className="ai-section">
      <div className="cli-head">
        <h3 className="section-title">{t('ai.title')}</h3>
        <span className="ai-launch">
          <button
            className="primary"
            disabled={!aiAvailable || poseCount === 0}
            title={aiAvailable ? undefined : availability?.reason}
            onClick={() => setDialogOpen(true)}
          >
            {t('ai.button')}
          </button>
          {availability && !aiAvailable && (
            <span className="muted ai-docs-hint" title={availability.reason}>
              {t('ai.unavailableDocs')}
            </span>
          )}
        </span>
      </div>
      {availability && !aiAvailable && availability.reason && (
        <p className="muted ai-reason">{availability.reason}</p>
      )}

      {runningTask && !isTerminalStatus(runningTask.status) && (
        <div className="ai-running">
          <div className="progress-head">
            <span className={`badge status-${runningTask.status}`}>
              {t(`status.${runningTask.status}`, runningTask.status)}
            </span>
            <span>{t('ai.running')}</span>
            <span className="muted">
              {runningTask.progress.percent}% · {runningTask.pose_indices.length}{' '}
              {t('ai.posesCount')}
            </span>
            <button className="danger-outline" onClick={() => void cancelTask(runningTask)}>
              {t('ai.cancel')}
            </button>
          </div>
          <div className="progress-bar">
            <div className="progress-fill" style={{ width: `${runningTask.progress.percent}%` }} />
          </div>
        </div>
      )}
      {runningTask && isTerminalStatus(runningTask.status) && runningTask.error && (
        <div className="banner">{runningTask.error}</div>
      )}
      {pollError && <p className="error-text">{pollError}</p>}

      {latestWithResults?.msa_server_used && (
        <p className="ai-warn">{t('ai.msaWarnHistory')}</p>
      )}
      {latestWithResults && <p className="muted">{t('ai.ligandScope')}</p>}

      {tasks.length > 0 && (
        <details className="ai-history" open={tasks.some((x) => !isTerminalStatus(x.status))}>
          <summary>
            {t('ai.history')}（{tasks.length}）
          </summary>
          <ul className="job-list">
            {[...tasks].reverse().map((task) => (
              <li key={task.id}>
                <div className="ai-task-row">
                  <span className={`badge status-${task.status}`}>
                    {t(`status.${task.status}`, task.status)}
                  </span>
                  <span className="muted">{new Date(task.created_at).toLocaleString()}</span>
                  <span className="muted">
                    {task.pose_indices.length} {t('ai.posesCount')}
                    {task.finished_at !== null && ` · ${t('status.completed')}`}
                  </span>
                  {!isTerminalStatus(task.status) && (
                    <button className="link-btn" onClick={() => void cancelTask(task)}>
                      {t('ai.cancel')}
                    </button>
                  )}
                </div>
              </li>
            ))}
          </ul>
        </details>
      )}

      {dialogOpen && (
        <AiRescoreDialog
          jobId={jobId}
          poseCount={poseCount}
          onClose={() => setDialogOpen(false)}
          onSubmitted={(task) => {
            setDialogOpen(false);
            setRunningTask(task);
            pushTasks([task, ...tasksRef.current]);
          }}
        />
      )}
    </div>
  );
}
