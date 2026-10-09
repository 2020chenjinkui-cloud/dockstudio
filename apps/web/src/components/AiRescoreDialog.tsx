import { useState } from 'react';
import { createAiRescore, getRescoreTask } from '../api';
import { useI18n } from '../i18n';
import { defaultPoseSelection, validateRescoreSelection } from '../rescore';
import type { AiRescoreTask } from '../types';

interface AiRescoreDialogProps {
  jobId: string;
  poseCount: number;
  onSubmitted: (task: AiRescoreTask) => void;
  onClose: () => void;
}

export default function AiRescoreDialog({
  jobId,
  poseCount,
  onSubmitted,
  onClose,
}: AiRescoreDialogProps) {
  const { t } = useI18n();
  const [selected, setSelected] = useState<number[]>(() => defaultPoseSelection(poseCount));
  const [useMsa, setUseMsa] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [busy, setBusy] = useState(false);

  const toggle = (poseIndex: number, on: boolean) => {
    setSelected((cur) =>
      on ? [...cur, poseIndex].sort((a, b) => a - b) : cur.filter((i) => i !== poseIndex),
    );
  };

  const submit = async () => {
    const validation = validateRescoreSelection(selected);
    if (!validation.ok) {
      setError(t('ai.noError'));
      return;
    }
    setBusy(true);
    setError(null);
    try {
      const res = await createAiRescore(jobId, {
        pose_indices: selected,
        use_msa_server: useMsa,
      });
      onSubmitted(await getRescoreTask(res.task_id));
    } catch (e) {
      setError(e instanceof Error ? e.message : String(e));
      setBusy(false);
    }
  };

  return (
    <div className="modal-overlay" onClick={onClose}>
      <div className="modal" onClick={(e) => e.stopPropagation()}>
        <h3>{t('ai.title')}</h3>
        <p className="muted">{t('ai.what')}</p>
        <p className="muted">{t('ai.poses')}</p>
        <div className="ai-pose-grid">
          {Array.from({ length: poseCount }, (_, i) => i + 1).map((poseIndex) => (
            <label key={poseIndex} className="check-row">
              <input
                type="checkbox"
                checked={selected.includes(poseIndex)}
                onChange={(e) => toggle(poseIndex, e.target.checked)}
              />
              #{poseIndex}
            </label>
          ))}
        </div>
        <label className="check-row">
          <input type="checkbox" checked={useMsa} onChange={(e) => setUseMsa(e.target.checked)} />
          {t('ai.msa')}
        </label>
        {useMsa && <p className="ai-warn">{t('ai.msaWarn')}</p>}
        {error && <p className="error-text">{error}</p>}
        <div className="modal-actions">
          <button onClick={onClose}>{t('bar.cancel')}</button>
          <button className="primary" onClick={() => void submit()} disabled={busy || selected.length === 0}>
            {busy ? t('ai.submitting') : t('ai.submit')}
          </button>
        </div>
      </div>
    </div>
  );
}
