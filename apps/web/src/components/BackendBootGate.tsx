import { useEffect, useRef, useState, type ReactNode } from 'react';
import { getHealth } from '../api';
import { useI18n } from '../i18n';

const POLL_MS = 2000;
const LONG_WAIT_S = 20;

interface BackendBootGateProps {
  children: ReactNode;
}

export default function BackendBootGate({ children }: BackendBootGateProps) {
  const { t } = useI18n();
  const [ready, setReady] = useState(false);
  const [elapsed, setElapsed] = useState(0);
  const [showLogs, setShowLogs] = useState(false);
  const timerRef = useRef<number | null>(null);

  useEffect(() => {
    if (ready) return;
    let cancelled = false;

    const attempt = () => {
      getHealth()
        .then(() => {
          if (!cancelled) setReady(true);
        })
        .catch(() => {
          if (!cancelled) {
            window.setTimeout(attempt, POLL_MS);
          }
        });
    };
    attempt();

    timerRef.current = window.setInterval(() => {
      setElapsed((s) => s + 1);
    }, 1000);

    return () => {
      cancelled = true;
      if (timerRef.current !== null) clearInterval(timerRef.current);
    };
  }, [ready]);

  if (ready) return <>{children}</>;

  const stage = elapsed < 6 ? 'gate.stage1' : elapsed < 12 ? 'gate.stage2' : 'gate.stage3';

  return (
    <div className="gate">
      <div className="card gate-card">
        <h2>{t('gate.title')}</h2>
        <div className="gate-spinner" />
        <p className="gate-stage">{t(stage)}</p>
        {elapsed >= LONG_WAIT_S && (
          <>
            <p className="muted">{t('gate.long')}</p>
            <div className="modal-actions gate-actions">
              <button onClick={() => setShowLogs((v) => !v)}>{t('gate.logs')}</button>
              <button
                className="primary"
                onClick={() => {
                  setElapsed(0);
                  setShowLogs(false);
                  getHealth()
                    .then(() => setReady(true))
                    .catch(() => undefined);
                }}
              >
                {t('gate.retry')}
              </button>
            </div>
            {showLogs && <p className="muted gate-log-text">{t('gate.logText')}</p>}
          </>
        )}
      </div>
    </div>
  );
}
