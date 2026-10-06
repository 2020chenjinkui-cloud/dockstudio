import { useEffect, useState } from 'react';
import { listPreparedLigands, prepareLigand } from '../../api';
import type { PreparedLigand } from '../../types';

interface StepLigandProps {
  ligand: PreparedLigand | null;
  onPrepared: (l: PreparedLigand) => void;
}

type Mode = 'file' | 'smiles' | 'library';

export default function StepLigand({ ligand, onPrepared }: StepLigandProps) {
  const [mode, setMode] = useState<Mode>('file');
  const [file, setFile] = useState<File | null>(null);
  const [smiles, setSmiles] = useState('');
  const [name, setName] = useState('');
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [library, setLibrary] = useState<PreparedLigand[] | null>(null);

  useEffect(() => {
    if (mode !== 'library' || library !== null) return;
    listPreparedLigands()
      .then(setLibrary)
      .catch((e: unknown) => setError(e instanceof Error ? e.message : String(e)));
  }, [mode, library]);

  const canSubmit = mode === 'file' ? file !== null : smiles.trim().length > 0;

  const submit = async () => {
    setBusy(true);
    setError(null);
    try {
      onPrepared(
        await prepareLigand({
          file: mode === 'file' ? (file ?? undefined) : undefined,
          smiles: mode === 'smiles' ? smiles.trim() : undefined,
          name: name.trim() || undefined,
        }),
      );
    } catch (e) {
      setError(e instanceof Error ? e.message : String(e));
    } finally {
      setBusy(false);
    }
  };

  return (
    <div className="step-body">
      <div className="tabs">
        <button className={mode === 'file' ? 'active' : ''} onClick={() => setMode('file')}>
          上传文件
        </button>
        <button className={mode === 'smiles' ? 'active' : ''} onClick={() => setMode('smiles')}>
          输入 SMILES
        </button>
        <button className={mode === 'library' ? 'active' : ''} onClick={() => setMode('library')}>
          已准备列表
        </button>
      </div>
      {mode !== 'library' && (
        <>
          {mode === 'file' ? (
            <input
              type="file"
              accept=".sdf,.mol,.mol2"
              onChange={(e) => setFile(e.target.files?.[0] ?? null)}
            />
          ) : (
            <input
              type="text"
              placeholder="例如：CC(=O)Oc1ccccc1C(=O)O"
              value={smiles}
              onChange={(e) => setSmiles(e.target.value)}
            />
          )}
          <div className="row">
            <input
              type="text"
              placeholder="名称（可选）"
              value={name}
              onChange={(e) => setName(e.target.value)}
            />
            <button onClick={submit} disabled={!canSubmit || busy}>
              {busy ? '处理中…' : '准备配体'}
            </button>
          </div>
        </>
      )}
      {mode === 'library' &&
        (!library ? (
          <p className="muted">加载中…</p>
        ) : library.length === 0 ? (
          <p className="muted">暂无已准备的配体。</p>
        ) : (
          <ul className="pick-list">
            {library.map((l) => (
              <li key={l.id}>
                <button
                  className={ligand?.id === l.id ? 'pick-item selected' : 'pick-item'}
                  onClick={() => onPrepared(l)}
                >
                  <span className="mono">{(l.name ?? l.id).slice(0, 24)}</span>
                  <span className="muted">可旋转键 {l.n_rotatable_bonds}</span>
                </button>
              </li>
            ))}
          </ul>
        ))}
      {error && <p className="error-text">{error}</p>}
      {ligand && (
        <dl className="kv">
          <dt>当前配体</dt>
          <dd className="mono">{ligand.name ?? ligand.id}</dd>
          <dt>来源</dt>
          <dd>{ligand.source === 'file' ? '文件' : 'SMILES'}</dd>
          <dt>可旋转键</dt>
          <dd>{ligand.n_rotatable_bonds}</dd>
          {ligand.smiles && (
            <>
              <dt>SMILES</dt>
              <dd className="mono wrap">{ligand.smiles}</dd>
            </>
          )}
        </dl>
      )}
    </div>
  );
}
