import { useEffect, useRef, useState } from 'react';
import { cancelBatch, createBatch, subscribeBatch } from '../../api';
import { useI18n } from '../../i18n';
import type { BatchSummary, EngineManifest } from '../../types';
import Tooltip from '../Tooltip';
import BatchPanel from './BatchPanel';
import type { WizardData } from './Wizard';

interface StepBatchProps {
  projectId: string;
  data: WizardData;
  engines: EngineManifest[];
  onBatchId: (batchId: string) => void;
  onOpenJob: (jobId: string) => void;
  onTerminal: (summary: BatchSummary) => void;
}

export default function StepBatch({
  projectId,
  data,
  engines,
  onBatchId,
  onOpenJob,
  onTerminal,
}: StepBatchProps) {
  const { t } = useI18n();
  const [submitting, setSubmitting] = useState(false);
  const [submitError, setSubmitError] = useState<string | null>(null);
  const [batchName, setBatchName] = useState('');
  const [summary, setSummary] = useState<BatchSummary | null>(null);
  const [pollError, setPollError] = useState<string | null>(null);
  const terminalRef = useRef(false);

  useEffect(() => {
    if (!data.batchId) return;
    setSummary(null);
    terminalRef.current = false;
    const stop = subscribeBatch(data.batchId, {
      onSummary: (s) => {
        setSummary(s);
        if (s.is_terminal && !terminalRef.current) {
          terminalRef.current = true;
          onTerminal(s);
        }
      },
      onError: (msg) => setPollError(msg),
    });
    return stop;
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [data.batchId]);

  const submit = async () => {
    if (!data.receptor || data.batchLigandIds.length === 0) return;
    setSubmitting(true);
    setSubmitError(null);
    try {
      const engine = engines.find((e) => e.available)?.name ?? 'vina';
      const res = await createBatch(projectId, {
        name: batchName.trim() || undefined,
        engine,
        receptor_id: data.receptor.id,
        ligand_ids: data.batchLigandIds,
        box: data.box,
        params: data.params,
      });
      onBatchId(res.batch_id);
    } catch (e) {
      setSubmitError(e instanceof Error ? e.message : String(e));
    } finally {
      setSubmitting(false);
    }
  };

  const cancel = async () => {
    if (!data.batchId) return;
    try {
      setSummary(await cancelBatch(data.batchId));
    } catch (e) {
      setPollError(e instanceof Error ? e.message : String(e));
    }
  };

  if (!data.batchId) {
    return (
      <div className="step-body">
        <dl className="kv">
          <dt>受体</dt>
          <dd className="mono">{data.receptor?.id.slice(0, 8) ?? '—'}</dd>
          <dt>配体数</dt>
          <dd>{data.batchLigandIds.length}</dd>
        </dl>
        <div className="control-row">
          <span>批量参数沿用上方参数卡（exhaustiveness {data.params.exhaustiveness} / num_modes {data.params.num_modes} / scoring {data.params.scoring}）</span>
          <Tooltip text={t('help.exhaustiveness')} />
        </div>
        <input
          type="text"
          placeholder="批次名称（可选）"
          value={batchName}
          onChange={(e) => setBatchName(e.target.value)}
        />
        {submitError && <div className="banner">{submitError}</div>}
        <button
          className="primary"
          onClick={submit}
          disabled={submitting || !data.receptor || data.batchLigandIds.length === 0}
        >
          {submitting ? '提交中…' : `提交批量任务（${data.batchLigandIds.length} 个配体）`}
        </button>
      </div>
    );
  }

  const seg = (n: number) =>
    summary && summary.total > 0 ? `${(n / summary.total) * 100}%` : '0%';

  return (
    <div className="step-body">
      <div className="progress-head">
        {summary ? (
          <>
            <span className={summary.is_terminal ? 'badge status-completed' : 'badge status-running'}>
              {summary.is_terminal ? '已完成' : '进行中'}
            </span>
            <span>
              完成 {summary.done}/{summary.total}
              {summary.failed > 0 && `，失败 ${summary.failed}`}
              {summary.canceled > 0 && `，取消 ${summary.canceled}`}
            </span>
            <span className="muted mono">batch {data.batchId.slice(0, 8)}</span>
          </>
        ) : (
          <span className="muted">连接中…</span>
        )}
      </div>
      {summary && (
        <div className="seg-bar" title={`done ${summary.done} / running ${summary.running} / failed ${summary.failed} / canceled ${summary.canceled} / queued ${summary.queued}`}>
          <span className="seg seg-done" style={{ width: seg(summary.done) }} />
          <span className="seg seg-running" style={{ width: seg(summary.running) }} />
          <span className="seg seg-failed" style={{ width: seg(summary.failed) }} />
          <span className="seg seg-canceled" style={{ width: seg(summary.canceled) }} />
          <span className="seg seg-queued" style={{ width: seg(summary.queued) }} />
        </div>
      )}
      {pollError && <p className="error-text">进度通道异常（已自动降级轮询）：{pollError}</p>}
      {summary && !summary.is_terminal && (
        <button className="danger-outline" onClick={cancel}>
          取消该批全部任务
        </button>
      )}
      {data.batchId && (
        <BatchPanel
          batchId={data.batchId}
          refreshKey={(summary?.done ?? 0) + (summary?.failed ?? 0) + (summary?.canceled ?? 0)}
          onOpenJob={onOpenJob}
        />
      )}
    </div>
  );
}
