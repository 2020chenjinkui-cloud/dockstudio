import { useEffect, useState } from 'react';
import { compareJobs } from '../api';
import { bestAffinity, toCompareCsv } from '../compare';
import { useI18n } from '../i18n';
import type { ComparisonRow, ComparisonResult } from '../types';

interface CompareDialogProps {
  jobIds: string[];
  onClose: () => void;
}

function shortSmiles(smiles: string | null, max = 40): string {
  if (!smiles) return '—';
  return smiles.length > max ? `${smiles.slice(0, max)}…` : smiles;
}

export default function CompareDialog({ jobIds, onClose }: CompareDialogProps) {
  const { t } = useI18n();
  const [result, setResult] = useState<ComparisonResult | null>(null);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    let cancelled = false;
    compareJobs(jobIds)
      .then((r) => {
        if (!cancelled) setResult(r);
      })
      .catch((e: unknown) => {
        if (!cancelled) setError(e instanceof Error ? e.message : String(e));
      });
    return () => {
      cancelled = true;
    };
  }, [jobIds]);

  const rows: ComparisonRow[] = result?.rows ?? [];
  const best = bestAffinity(rows);

  const exportCsv = () => {
    const csv = toCompareCsv(rows);
    const blob = new Blob([`﻿${csv}`], { type: 'text/csv;charset=utf-8' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = 'compare.csv';
    a.click();
    URL.revokeObjectURL(url);
  };

  return (
    <div className="modal-overlay" onClick={onClose}>
      <div className="modal compare-modal" onClick={(e) => e.stopPropagation()}>
        <div className="cli-head">
          <h3>
            {t('cmp.title')}（{jobIds.length}）
          </h3>
          <button onClick={onClose}>{t('bar.close')}</button>
        </div>
        {error && <div className="banner">{error}</div>}
        {!result && !error && <p className="muted">{t('common.loading')}</p>}
        {result && result.rows.length === 0 && (
          <p className="muted">{t('cmp.hint')}</p>
        )}
        {result && result.rows.length > 0 && (
          <>
            {result.rows.length !== jobIds.length && (
              <p className="muted">
                {t('cmp.hint')}（{result.rows.length}/{jobIds.length}）
              </p>
            )}
            <div className="table-scroll">
              <table className="data-table">
                <thead>
                  <tr>
                    <th>{t('cmp.col.ligand')}</th>
                    <th>{t('cmp.col.smiles')}</th>
                    <th>{t('cmp.col.engine')}</th>
                    <th>{t('cmp.col.scoring')}</th>
                    <th>{t('cmp.col.exhaustiveness')}</th>
                    <th>{t('cmp.col.affinity')}</th>
                    <th>{t('cmp.col.poses')}</th>
                    <th>{t('cmp.col.status')}</th>
                  </tr>
                </thead>
                <tbody>
                  {rows.map((r) => (
                    <tr key={r.job_id}>
                      <td className="mono">{r.ligand_name || r.ligand_id.slice(0, 8)}</td>
                      <td className="mono wrap smiles-cell">{shortSmiles(r.smiles)}</td>
                      <td>{r.engine}</td>
                      <td>{r.scoring}</td>
                      <td>{r.exhaustiveness}</td>
                      <td className={best !== null && r.top_affinity === best ? 'best-cell' : ''}>
                        {r.top_affinity !== null ? r.top_affinity.toFixed(2) : '—'}
                        {best !== null && r.top_affinity === best && (
                          <span className="muted">（{t('cmp.best')}）</span>
                        )}
                      </td>
                      <td>{r.n_poses}</td>
                      <td>
                        <span className={`badge status-${r.status}`}>
                          {t(`status.${r.status}`, r.status)}
                        </span>
                      </td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
            <div className="modal-actions">
              <button className="primary" onClick={exportCsv}>
                {t('cmp.export')}
              </button>
            </div>
          </>
        )}
      </div>
    </div>
  );
}
