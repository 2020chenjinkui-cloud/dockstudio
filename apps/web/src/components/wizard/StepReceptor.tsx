import { useEffect, useState } from 'react';
import { listPreparedReceptors, prepareReceptor } from '../../api';
import type { PreparedReceptor } from '../../types';
import FormatHint from '../FormatHint';

interface StepReceptorProps {
  receptor: PreparedReceptor | null;
  onPrepared: (r: PreparedReceptor) => void;
  onAutoBox: (center: [number, number, number]) => void;
}

type Mode = 'upload' | 'library';

export default function StepReceptor({ receptor, onPrepared, onAutoBox }: StepReceptorProps) {
  const [mode, setMode] = useState<Mode>('upload');
  const [file, setFile] = useState<File | null>(null);
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [library, setLibrary] = useState<PreparedReceptor[] | null>(null);
  const [selectedCocrystal, setSelectedCocrystal] = useState<number | null>(null);

  useEffect(() => {
    if (mode !== 'library' || library !== null) return;
    listPreparedReceptors()
      .then(setLibrary)
      .catch((e: unknown) => setError(e instanceof Error ? e.message : String(e)));
  }, [mode, library]);

  const upload = async () => {
    if (!file) return;
    setBusy(true);
    setError(null);
    setSelectedCocrystal(null);
    try {
      onPrepared(await prepareReceptor(file));
    } catch (e) {
      setError(e instanceof Error ? e.message : String(e));
    } finally {
      setBusy(false);
    }
  };

  return (
    <div className="step-body">
      <FormatHint kind="receptor" />
      <div className="tabs">
        <button className={mode === 'upload' ? 'active' : ''} onClick={() => setMode('upload')}>
          上传 PDB
        </button>
        <button className={mode === 'library' ? 'active' : ''} onClick={() => setMode('library')}>
          已准备列表
        </button>
      </div>
      {mode === 'upload' ? (
        <div className="row">
          <input
            type="file"
            accept=".pdb,.cif,.mmcif"
            onChange={(e) => setFile(e.target.files?.[0] ?? null)}
          />
          <button onClick={upload} disabled={!file || busy}>
            {busy ? '处理中…' : '上传并准备'}
          </button>
        </div>
      ) : !library ? (
        <p className="muted">加载中…</p>
      ) : library.length === 0 ? (
        <p className="muted">暂无已准备的受体，请先上传。</p>
      ) : (
        <ul className="pick-list">
          {library.map((r) => (
            <li key={r.id}>
              <button
                className={receptor?.id === r.id ? 'pick-item selected' : 'pick-item'}
                onClick={() => onPrepared(r)}
              >
                <span className="mono">{r.id.slice(0, 8)}</span>
                <span className="muted">链 {r.chains.join(', ') || '—'}</span>
                {r.parent_id != null ? (
                  <span className="badge badge-flex">柔性：{r.flex_residues.join(' ')}</span>
                ) : (
                  <span className="muted">共晶配体 {r.cocrystal_ligands.length} 个</span>
                )}
              </button>
            </li>
          ))}
        </ul>
      )}
      {error && <p className="error-text">{error}</p>}
      {receptor && (
        <div className="summary">
          <dl className="kv">
            <dt>当前受体</dt>
            <dd className="mono">{receptor.id}</dd>
            <dt>链</dt>
            <dd>{receptor.chains.length > 0 ? receptor.chains.join(', ') : '（无）'}</dd>
            {receptor.parent_id != null && (
              <>
                <dt>柔性残基</dt>
                <dd>
                  <span className="badge badge-flex">{receptor.flex_residues.join(' ')}</span>
                  <span className="muted">（基于 {receptor.parent_id.slice(0, 8)}）</span>
                </dd>
              </>
            )}
          </dl>
          {receptor.warnings.length > 0 && (
            <ul className="warnings">
              {receptor.warnings.map((w, i) => (
                <li key={i}>{w}</li>
              ))}
            </ul>
          )}
          {receptor.cocrystal_ligands.length > 0 ? (
            <table className="data-table">
              <thead>
                <tr>
                  <th>共晶配体</th>
                  <th>链</th>
                  <th>原子数</th>
                  <th></th>
                </tr>
              </thead>
              <tbody>
                {receptor.cocrystal_ligands.map((lig, i) => (
                  <tr key={i} className={selectedCocrystal === i ? 'selected' : ''}>
                    <td className="mono">{lig.resname}</td>
                    <td>{lig.chain}</td>
                    <td>{lig.n_atoms}</td>
                    <td>
                      <button
                        className="link-btn"
                        onClick={() => {
                          setSelectedCocrystal(i);
                          onAutoBox(lig.center);
                        }}
                      >
                        定盒
                      </button>
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          ) : (
            <p className="muted">未检测到共晶配体，下一步需手动填写盒子中心。</p>
          )}
        </div>
      )}
    </div>
  );
}
