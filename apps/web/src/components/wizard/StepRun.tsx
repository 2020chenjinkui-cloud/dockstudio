import { useEffect, useMemo, useRef, useState } from 'react';
import { cancelJob, createFlexReceptor, createJob, subscribeJob } from '../../api';
import { parseFlexInput, validateFlexResidues } from '../../flex';
import { PARAM_HELP } from '../../helpText';
import { isTerminalStatus, stageLabel, statusLabel } from '../../labels';
import type { DockParams, EngineManifest, Job } from '../../types';
import Tooltip from '../Tooltip';
import type { WizardData } from './Wizard';

interface StepRunProps {
  projectId: string;
  data: WizardData;
  advanced: boolean;
  engines: EngineManifest[];
  onParamsChange: (params: DockParams) => void;
  onFlexText: (text: string) => void;
  onJobId: (jobId: string) => void;
  onJobUpdate: (job: Job) => void;
  onCompleted: () => void;
}

export default function StepRun({
  projectId,
  data,
  advanced,
  engines,
  onParamsChange,
  onFlexText,
  onJobId,
  onJobUpdate,
  onCompleted,
}: StepRunProps) {
  const [submitting, setSubmitting] = useState(false);
  const [submitError, setSubmitError] = useState<string | null>(null);
  const [flexNote, setFlexNote] = useState<string | null>(null);
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

  const flexValidation = useMemo(() => {
    if (!advanced || data.flexText.trim() === '') return { valid: [], errors: [] };
    return validateFlexResidues(parseFlexInput(data.flexText), data.receptor?.chains ?? []);
  }, [advanced, data.flexText, data.receptor]);

  const flexBlocked =
    advanced && data.flexText.trim() !== '' && flexValidation.errors.length > 0;

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

  const num = (v: string): number => {
    const n = Number.parseFloat(v);
    return Number.isNaN(n) ? 0 : n;
  };

  const submit = async () => {
    if (!data.receptor || !data.ligand || flexBlocked) return;
    setSubmitting(true);
    setSubmitError(null);
    setFlexNote(null);
    try {
      const engine = engines.find((e) => e.available)?.name ?? 'vina';
      let receptorId = data.receptor.id;
      if (advanced && flexValidation.valid.length > 0) {
        const variant = await createFlexReceptor(data.receptor.id, flexValidation.valid);
        receptorId = variant.id;
        setFlexNote(
          `已创建柔性变体 ${variant.id.slice(0, 8)}（${flexValidation.valid.join(' ')}）`,
        );
      }
      const res = await createJob(projectId, {
        engine,
        receptor_id: receptorId,
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

  if (!data.jobId) {
    return (
      <div className="step-body">
        <div className="grid-form">
          <label>
            <span className="control-row">
              exhaustiveness (1-64)
              <Tooltip text={PARAM_HELP.exhaustiveness} />
            </span>
            <input
              type="number"
              min={1}
              max={64}
              value={data.params.exhaustiveness}
              onChange={(e) =>
                onParamsChange({
                  ...data.params,
                  exhaustiveness: Math.min(64, Math.max(1, num(e.target.value))),
                })
              }
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
          {advanced && (
            <>
              <label>
                <span className="control-row">
                  energy_range
                  <Tooltip text={PARAM_HELP.energyRange} />
                </span>
                <input
                  type="number"
                  step={0.5}
                  min={1}
                  value={data.params.energy_range}
                  onChange={(e) =>
                    onParamsChange({ ...data.params, energy_range: num(e.target.value) })
                  }
                />
              </label>
              <label>
                <span className="control-row">
                  seed（留空随机）
                  <Tooltip text={PARAM_HELP.seed} />
                </span>
                <input
                  type="text"
                  placeholder="随机"
                  value={data.params.seed === null ? '' : String(data.params.seed)}
                  onChange={(e) => {
                    const t = e.target.value.trim();
                    onParamsChange({
                      ...data.params,
                      seed: t === '' ? null : Math.trunc(num(t)),
                    });
                  }}
                />
              </label>
            </>
          )}
        </div>
        {advanced && (
          <div className="flex-card">
            <div className="control-row">
              <span>柔性残基（最多 10 个）</span>
              <Tooltip text={PARAM_HELP.flexResidues} />
            </div>
            <input
              type="text"
              placeholder="如：A:42 A:57 B:101（留空为刚性对接）"
              value={data.flexText}
              onChange={(e) => onFlexText(e.target.value)}
            />
            {flexValidation.errors.length > 0 && (
              <ul className="error-list">
                {flexValidation.errors.map((err, i) => (
                  <li key={i}>{err}</li>
                ))}
              </ul>
            )}
            {flexValidation.errors.length === 0 && flexValidation.valid.length > 0 && (
              <p className="muted">
                将对 {flexValidation.valid.length} 个残基启用柔性：{flexValidation.valid.join(' ')}
              </p>
            )}
          </div>
        )}
        <p className="muted">预计等待数十秒到数分钟，取决于 exhaustiveness、盒子大小与本机性能。</p>
        {flexNote && <p className="ok-text">{flexNote}</p>}
        {submitError && <div className="banner">{submitError}</div>}
        <button
          className="primary"
          onClick={submit}
          disabled={submitting || !data.receptor || !data.ligand || flexBlocked}
        >
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
      {flexNote && <p className="ok-text">{flexNote}</p>}
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
