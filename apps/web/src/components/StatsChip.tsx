import { useEffect, useState } from 'react';
import { getStats } from '../api';
import { useI18n } from '../i18n';
import type { Stats } from '../types';

export default function StatsChip() {
  const { t } = useI18n();
  const [stats, setStats] = useState<Stats | null>(null);

  useEffect(() => {
    let cancelled = false;
    getStats()
      .then((s) => {
        if (!cancelled) setStats(s);
      })
      .catch(() => {
        // 后端无 /api/stats（旧版）时静默隐藏
      });
    return () => {
      cancelled = true;
    };
  }, []);

  if (!stats || (stats.projects === 0 && stats.completed_jobs === 0)) return null;

  return (
    <span className="stats-chip" title={Object.entries(stats.by_engine).map(([k, v]) => `${k}: ${v}`).join(' / ')}>
      {t('stats.projects')} {stats.projects} · {t('stats.completed')} {stats.completed_jobs}
      {stats.avg_dock_seconds !== null &&
        ` · ${t('stats.avg')} ${stats.avg_dock_seconds.toFixed(1)} s`}
    </span>
  );
}
