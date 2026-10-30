import { useEffect, useMemo, useRef, useState } from 'react';
import { cancelJob, createFlexReceptor, createJob, subscribeJob } from '../../api';
import { parseFlexInput, validateFlexResidues } from '../../flex';
import { useI18n } from '../../i18n';
import { isTerminalStatus } from '../../labels';
import type { DockParams, EngineManifest, Job } from '../../types';
import Tooltip from '../Tooltip';
import ParamsPanel from './ParamsPanel';
import type { WizardData } from './Wizard';

interface StepRunProps {
  projectId: string;
  data: WizardData;
  advanced: boolean;
  engines: EngineManifest[];
  onParamsChange: (params: DockParams) => void;
  onEngineChange: (engine: string) => void;
  onFlexText: (text: string) => void;
  onJobId: (jobId: string) => void;
  onJobUpdate: (job: Job) => void;
  onCompleted: () => void;
  onReset: () => void;
}

export default function StepRun({
  projectId,
  data,
  advanced,
  engines,
  onParamsChange,
  onEngineChange,
  onFlexText,
  onJobId,
  onJobUpdate,
  onCompleted,
  onReset,
}: StepRunProps) {
  const { t } = useI18n();
  const [submitting, setSubmitting] = useState(false);
  const [submitError, setSubmitError] = useState<string | null>(null);
  const [flexNote, setFlexNote] = useState<string | null>(null);
  const [job, setJob] = useState<Job | null>(null);
  const [pollError, setPollError] = useState<string | null>(null);
  const [beats, setBeats] = useState(0);
  const completedRef = useRef(false);

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

  const submit = async () => {
    if (!data.receptor || !data.ligand || flexBlocked) return;
    setSubmitting(true);
    setSubmitError(null);
    setFlexNote(null);
    try {
      let receptorId = data.receptor.id;
      if (advanced && flexValidation.valid.length > 0) {
        const variant = await createFlexReceptor(data.receptor.id, flexValidation.valid);
        receptorId = variant.id;
        setFlexNote(
          `已创建柔性变体 ${variant.id.slice(0, 8)}（${flexValidation.valid.join(' ')}）`,
        );
      }
      const res = await createJob(projectId, {
        engine: data.engine,
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
        <ParamsPanel
          params={data.params}
          engines={engines}
          engine={data.engine}
          advanced={advanced}
          onEngineChange={onEngineChange}
          onParamsChange={onParamsChange}
        />
        {advanced && (
          <div className="flex-card">
            <div className="control-row">
              <span>柔性残基（最多 10 个）</span>
              <Tooltip text={t('help.flexResidues')} />
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
        <span className={`badge status-${status}`}>{t(`status.${status}`, status)}</span>
        <span>{progress ? t(`stage.${progress.stage}`, progress.stage) : '连接中…'}</span>
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
      {terminal && (
        <button className="primary" onClick={onReset}>
          {t('nav.newRun')}
        </button>
      )}
    </div>
  );
}
