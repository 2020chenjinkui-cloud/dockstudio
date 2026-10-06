import { useState } from 'react';
import { prepareLigand } from '../api';
import type { PreparedLigand } from '../types';

interface LigandCardProps {
  ligand: PreparedLigand | null;
  onPrepared: (ligand: PreparedLigand) => void;
}

type Tab = 'file' | 'smiles';

export default function LigandCard({ ligand, onPrepared }: LigandCardProps) {
  const [tab, setTab] = useState<Tab>('file');
  const [file, setFile] = useState<File | null>(null);
  const [smiles, setSmiles] = useState('');
  const [name, setName] = useState('');
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const canSubmit = tab === 'file' ? file !== null : smiles.trim().length > 0;

  const upload = async () => {
    setBusy(true);
    setError(null);
    try {
      onPrepared(
        await prepareLigand({
          file: tab === 'file' ? (file ?? undefined) : undefined,
          smiles: tab === 'smiles' ? smiles.trim() : undefined,
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
    <section className="card">
      <h2>2. 配体</h2>
      <div className="tabs">
        <button className={tab === 'file' ? 'active' : ''} onClick={() => setTab('file')}>
          上传文件
        </button>
        <button className={tab === 'smiles' ? 'active' : ''} onClick={() => setTab('smiles')}>
          输入 SMILES
        </button>
      </div>
      {tab === 'file' ? (
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
        <button onClick={upload} disabled={!canSubmit || busy}>
          {busy ? '处理中…' : '准备配体'}
        </button>
      </div>
      {error && <p className="error-text">{error}</p>}
      {ligand && (
        <dl className="kv">
          <dt>ID</dt>
          <dd className="mono">{ligand.id}</dd>
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
    </section>
  );
}
