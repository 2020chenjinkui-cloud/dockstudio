import { useCallback, useEffect, useState } from 'react';
import { batchResultsCsvUrl, getBatchResults } from '../../api';
import { useI18n } from '../../i18n';
import type { BatchResults } from '../../types';

const PAGE_SIZE = 20;

interface BatchPanelProps {
  batchId: string;
  refreshKey?: number;
  onOpenJob: (jobId: string) => void;
}

export default function BatchPanel({ batchId, refreshKey = 0, onOpenJob }: BatchPanelProps) {
  const { t } = useI18n();
  const [results, setResults] = useState<BatchResults | null>(null);
  const [offset, setOffset] = useState(0);
  const [error, setError] = useState<string | null>(null);

  const load = useCallback(async () => {
    try {
      const data = await getBatchResults(batchId, {
        sort: 'affinity',
        order: 'asc',
        limit: PAGE_SIZE,
        offset,
      });
      setResults(data);
      setError(null);
    } catch (e) {
      setError(e instanceof Error ? e.message : String(e));
    }
  }, [batchId, offset]);

  useEffect(() => {
    void load();
  }, [load, refreshKey]);

  const total = results?.total ?? 0;
  const pageCount = Math.max(1, Math.ceil(total / PAGE_SIZE));
  const page = Math.floor(offset / PAGE_SIZE) + 1;

  return (
    <div className="batch-panel">
      <div className="batch-toolbar">
        <h3 className="section-title">排行榜（按 top affinity 排序）</h3>
        <span className="batch-toolbar-actions">
          <button onClick={() => void load()}>刷新</button>
          <a className="btn-link small" href={batchResultsCsvUrl(batchId)} download>
            导出 CSV
          </a>
        </span>
      </div>
      {error && <p className="error-text">{error}</p>}
      {!results ? (
        <p className="muted">加载中…</p>
      ) : results.rows.length === 0 ? (
        <p className="muted">尚无已完成的条目。</p>
      ) : (
        <table className="data-table">
          <thead>
            <tr>
              <th>#</th>
              <th>配体</th>
              <th>SMILES</th>
              <th>状态</th>
              <th>top affinity</th>
              <th>错误</th>
            </tr>
          </thead>
          <tbody>
            {results.rows.map((row, i) => (
              <tr key={row.job_id} onClick={() => onOpenJob(row.job_id)} title="查看该任务结果">
                <td>{offset + i + 1}</td>
                <td className="mono">{row.ligand_name || row.ligand_id.slice(0, 8)}</td>
                <td className="mono wrap smiles-cell">{row.smiles}</td>
                <td>
                  <span className={`badge status-${row.status}`}>{t(`status.${row.status}`, row.status)}</span>
                </td>
                <td>{row.top_affinity !== null ? row.top_affinity.toFixed(1) : '—'}</td>
                <td className="error-text cell-error">{row.error ?? ''}</td>
              </tr>
            ))}
          </tbody>
        </table>
      )}
      {total > PAGE_SIZE && (
        <div className="pagination">
          <button disabled={offset === 0} onClick={() => setOffset(Math.max(0, offset - PAGE_SIZE))}>
            上一页
          </button>
          <span className="muted">
            第 {page} / {pageCount} 页（共 {total} 条）
          </span>
          <button
            disabled={offset + PAGE_SIZE >= total}
            onClick={() => setOffset(offset + PAGE_SIZE)}
          >
            下一页
          </button>
        </div>
      )}
    </div>
  );
}
