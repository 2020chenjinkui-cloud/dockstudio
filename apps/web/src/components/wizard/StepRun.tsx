import { useEffect, useMemo, useRef, useState } from 'react';
import { cancelJob, createJob, subscribeJob } from '../../api';
import { PARAM_HELP } from '../../helpText';
import { isTerminalStatus, stageLabel, statusLabel } from '../../labels';
import type { DockParams, EngineManifest, Job } from '../../types';
import Tooltip from '../Tooltip';
import type { WizardData } from './Wizard';

interface StepRunProps {
  projectId: string;
  data: WizardData;
  engines: EngineManifest[];
  onParamsChange: (params: DockParams) => void;
  onJobId: (jobId: string) => void;
  onJobUpdate: (job: Job) => void;
  onCompleted: () => void;
}

export default function StepRun({
  projectId,
  data,
  engines,
  onParamsChange,
  onJobId,
  onJobUpdate,
  onCompleted,
}: StepRunProps) {
  const [submitting, setSubmitting] = useState(false);
  const [submitError, setSubmitError] = useState<string | null>(null);
  const [job, setJob] = useState<Job | null>(null);
  const [pollError, setPollError] = useState<string | null>(null);
  const [beats, setBeats] = useState(0);
  const completedRef = useRef(false);

  const scoringOptions = useMemo(
    () => [...new Set(engines.flatMap((e) => e.scoring_functions))],
    [engines],
  );
  const scoringList = scoringOptions.includes(data.params.scoring)
    ? scoringOptions
    : [data.params.scoring, ...scoringOptions];

  useEffect(() => {
    if (!data.jobId) return;
    setJob(null);
    setBeats(0);
    completedRef.current = false;
    const stop = subscribeJob(data.jobId, {
      onProgress: () => setBeats((n) => n + 1),
      onJob: (j) => {
        setJob(j);
        onJobUpdate(j);
        if (j.status === 'completed' && !completedRef.current) {
          completedRef.current = true;
          onCompleted();
        }
      },
      onError: (msg) => setPollError(msg),
    });
    return stop;
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [data.jobId]);

  const submit = async () => {
    if (!data.receptor || !data.ligand) return;
    setSubmitting(true);
    setSubmitError(null);
    try {
      const engine = engines.find((e) => e.available)?.name ?? 'vina';
      const res = await createJob(projectId, {
        engine,
        receptor_id: data.receptor.id,
        ligand_id: data.ligand.id,
        box: data.box,
        params: data.params,
      });
      onJobId(res.job_id);
    } catch (e) {
      setSubmitError(e instanceof Error ? e.message : String(e));
    } finally {
      setSubmitting(false);
    }
  };

  const cancel = async () => {
    if (!data.jobId) return;
    try {
      const j = await cancelJob(data.jobId);
      setJob(j);
      onJobUpdate(j);
    } catch (e) {
      setPollError(e instanceof Error ? e.message : String(e));
    }
  };

  const num = (v: string): number => {
    const n = Number.parseFloat(v);
    return Number.isNaN(n) ? 0 : n;
  };

  if (!data.jobId) {
    return (
      <div className="step-body">
        <div className="grid-form">
          <label>
            <span className="control-row">
              exhaustiveness
              <Tooltip text={PARAM_HELP.exhaustiveness} />
            </span>
            <input
              type="number"
              min={1}
              value={data.params.exhaustiveness}
              onChange={(e) => onParamsChange({ ...data.params, exhaustiveness: num(e.target.value) })}
            />
          </label>
          <label>
            <span className="control-row">
              num_modes
              <Tooltip text={PARAM_HELP.num_modes} />
            </span>
            <input
              type="number"
              min={1}
              value={data.params.num_modes}
              onChange={(e) => onParamsChange({ ...data.params, num_modes: num(e.target.value) })}
            />
          </label>
          <label>
            <span className="control-row">
              cpu
              <Tooltip text={PARAM_HELP.cpu} />
            </span>
            <input
              type="number"
              min={1}
              value={data.params.cpu}
              onChange={(e) => onParamsChange({ ...data.params, cpu: num(e.target.value) })}
            />
          </label>
          <label>
            <span className="control-row">
              scoring
              <Tooltip text={PARAM_HELP.scoring} />
            </span>
            <select
              value={data.params.scoring}
              onChange={(e) => onParamsChange({ ...data.params, scoring: e.target.value })}
            >
              {scoringList.map((s) => (
                <option key={s} value={s}>
                  {s}
                </option>
              ))}
            </select>
          </label>
        </div>
        <p className="muted">预计等待数十秒到数分钟，取决于 exhaustiveness、盒子大小与本机性能。</p>
        {submitError && <div className="banner">{submitError}</div>}
        <button className="primary" onClick={submit} disabled={submitting || !data.receptor || !data.ligand}>
          {submitting ? '提交中…' : '提交对接任务'}
        </button>
      </div>
    );
  }

  const status = job?.status ?? 'queued';
  const progress = job?.progress;
  const terminal = isTerminalStatus(status);

  return (
    <div className="step-body">
      <div className="progress-head">
        <span className={`badge status-${status}`}>{statusLabel(status)}</span>
        <span>{progress ? stageLabel(progress.stage) : '连接中…'}</span>
        <span className="muted mono">job {data.jobId.slice(0, 8)}</span>
      </div>
      <div className="progress-bar">
        <div className="progress-fill" style={{ width: `${progress?.percent ?? 0}%` }} />
      </div>
      <p className="muted">
        {progress?.message ?? '等待后端推送进度…'}（已接收 {beats} 条进度）
      </p>
      {pollError && <p className="error-text">进度通道异常（已自动降级轮询）：{pollError}</p>}
      {job?.error && <div className="banner">{job.error}</div>}
      {status === 'failed' && (
        <p className="muted">任务失败，可返回前面的步骤修改参数后重新提交。</p>
      )}
      {!terminal && (
        <button className="danger-outline" onClick={cancel}>
          取消任务
        </button>
      )}
    </div>
  );
}
