import { useState } from 'react';
import { prepareReceptor } from '../api';
import type { PreparedReceptor } from '../types';

interface ReceptorCardProps {
  receptor: PreparedReceptor | null;
  onPrepared: (receptor: PreparedReceptor) => void;
  onSelectCocrystal: (center: [number, number, number]) => void;
}

export default function ReceptorCard({ receptor, onPrepared, onSelectCocrystal }: ReceptorCardProps) {
  const [file, setFile] = useState<File | null>(null);
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [selectedIdx, setSelectedIdx] = useState<number | null>(null);

  const upload = async () => {
    if (!file) return;
    setBusy(true);
    setError(null);
    setSelectedIdx(null);
    try {
      onPrepared(await prepareReceptor(file));
    } catch (e) {
      setError(e instanceof Error ? e.message : String(e));
    } finally {
      setBusy(false);
    }
  };

  return (
    <section className="card">
      <h2>1. 受体</h2>
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
      {error && <p className="error-text">{error}</p>}
      {receptor && (
        <>
          <dl className="kv">
            <dt>ID</dt>
            <dd className="mono">{receptor.id}</dd>
            <dt>链</dt>
            <dd>{receptor.chains.length > 0 ? receptor.chains.join(', ') : '（无）'}</dd>
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
                  <tr key={i} className={selectedIdx === i ? 'selected' : ''}>
                    <td className="mono">{lig.resname}</td>
                    <td>{lig.chain}</td>
                    <td>{lig.n_atoms}</td>
                    <td>
                      <button
                        className="link-btn"
                        onClick={() => {
                          setSelectedIdx(i);
                          onSelectCocrystal(lig.center);
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
            <p className="muted">未检测到共晶配体，请手动填写盒子中心。</p>
          )}
        </>
      )}
    </section>
  );
}
