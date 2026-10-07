import { useState } from 'react';
import { reportUrlV2 } from '../api';
import { useI18n } from '../i18n';
import { REPORT_SECTIONS } from '../types';

interface ReportDialogProps {
  jobId: string;
  onClose: () => void;
}

export default function ReportDialog({ jobId, onClose }: ReportDialogProps) {
  const { t } = useI18n();
  const [checked, setChecked] = useState<string[]>([...REPORT_SECTIONS]);
  const [title, setTitle] = useState('');

  const toggle = (section: string, on: boolean) => {
    setChecked((cur) => (on ? [...cur, section] : cur.filter((s) => s !== section)));
  };

  const download = () => {
    window.location.assign(reportUrlV2(jobId, { sections: checked, title }));
    onClose();
  };

  return (
    <div className="modal-overlay" onClick={onClose}>
      <div className="modal" onClick={(e) => e.stopPropagation()}>
        <h3>{t('report.dialogTitle')}</h3>
        <p className="muted">{t('report.sections')}</p>
        <div className="section-checks">
          {REPORT_SECTIONS.map((s) => (
            <label key={s} className="check-row">
              <input
                type="checkbox"
                checked={checked.includes(s)}
                onChange={(e) => toggle(s, e.target.checked)}
              />
              {t(`report.section.${s}`)}
            </label>
          ))}
        </div>
        <input
          type="text"
          maxLength={80}
          placeholder={t('report.customTitle')}
          value={title}
          onChange={(e) => setTitle(e.target.value)}
        />
        <div className="modal-actions">
          <button onClick={onClose}>{t('bar.cancel')}</button>
          <button className="primary" onClick={download} disabled={checked.length === 0}>
            {t('report.confirm')}
          </button>
        </div>
      </div>
    </div>
  );
}
