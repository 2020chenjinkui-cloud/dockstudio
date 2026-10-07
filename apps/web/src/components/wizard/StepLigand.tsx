import { useEffect, useState } from 'react';
import { listPreparedLigands, prepareLibrary, prepareLigand } from '../../api';
import { PARAM_HELP } from '../../helpText';
import type { PreparedLibrary, PreparedLigand } from '../../types';
import Tooltip from '../Tooltip';

interface StepLigandProps {
  mode: 'single' | 'batch';
  ligand: PreparedLigand | null;
  onPrepared: (l: PreparedLigand) => void;
  batchLigandIds: string[];
  onBatchChange: (ids: string[]) => void;
}

type Mode = 'file' | 'smiles' | 'library';
type BatchModeTab = 'library-file' | 'multi-select';

export default function StepLigand({
  mode,
  ligand,
  onPrepared,
  batchLigandIds,
  onBatchChange,
}: StepLigandProps) {
  if (mode === 'batch') {
    return <BatchLigands batchLigandIds={batchLigandIds} onBatchChange={onBatchChange} />;
  }
  return (
    <SingleLigand ligand={ligand} onPrepared={onPrepared} />
  );
}

function SingleLigand({
  ligand,
  onPrepared,
}: {
  ligand: PreparedLigand | null;
  onPrepared: (l: PreparedLigand) => void;
}) {
  const [tab, setTab] = useState<Mode>('file');
  const [file, setFile] = useState<File | null>(null);
  const [smiles, setSmiles] = useState('');
  const [name, setName] = useState('');
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [library, setLibrary] = useState<PreparedLigand[] | null>(null);

  useEffect(() => {
    if (tab !== 'library' || library !== null) return;
    listPreparedLigands()
      .then(setLibrary)
      .catch((e: unknown) => setError(e instanceof Error ? e.message : String(e)));
  }, [tab, library]);

  const canSubmit = tab === 'file' ? file !== null : smiles.trim().length > 0;

  const submit = async () => {
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
    <div className="step-body">
      <div className="tabs">
        <button className={tab === 'file' ? 'active' : ''} onClick={() => setTab('file')}>
          上传文件
        </button>
        <button className={tab === 'smiles' ? 'active' : ''} onClick={() => setTab('smiles')}>
          输入 SMILES
        </button>
        <button className={tab === 'library' ? 'active' : ''} onClick={() => setTab('library')}>
          已准备列表
        </button>
      </div>
      {tab !== 'library' && (
        <>
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
            <button onClick={submit} disabled={!canSubmit || busy}>
              {busy ? '处理中…' : '准备配体'}
            </button>
          </div>
        </>
      )}
      {tab === 'library' &&
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

function BatchLigands({
  batchLigandIds,
  onBatchChange,
}: {
  batchLigandIds: string[];
  onBatchChange: (ids: string[]) => void;
}) {
  const [tab, setTab] = useState<BatchModeTab>('library-file');
  const [file, setFile] = useState<File | null>(null);
  const [libName, setLibName] = useState('');
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [preparedLib, setPreparedLib] = useState<PreparedLibrary | null>(null);
  const [library, setLibrary] = useState<PreparedLigand[] | null>(null);

  useEffect(() => {
    if (tab !== 'multi-select' || library !== null) return;
    listPreparedLigands()
      .then(setLibrary)
      .catch((e: unknown) => setError(e instanceof Error ? e.message : String(e)));
  }, [tab, library]);

  const upload = async () => {
    if (!file) return;
    setBusy(true);
    setError(null);
    setPreparedLib(null);
    try {
      const lib = await prepareLibrary(file, libName.trim() || undefined);
      setPreparedLib(lib);
      onBatchChange(lib.ligand_ids);
    } catch (e) {
      setError(e instanceof Error ? e.message : String(e));
    } finally {
      setBusy(false);
    }
  };

  const toggle = (id: string, checked: boolean) => {
    onBatchChange(
      checked ? [...batchLigandIds, id] : batchLigandIds.filter((x) => x !== id),
    );
  };

  return (
    <div className="step-body">
      <div className="tabs">
        <button
          className={tab === 'library-file' ? 'active' : ''}
          onClick={() => setTab('library-file')}
        >
          上传库文件
        </button>
        <button
          className={tab === 'multi-select' ? 'active' : ''}
          onClick={() => setTab('multi-select')}
        >
          从已准备配体多选
        </button>
      </div>
      {tab === 'library-file' ? (
        <>
          <div className="control-row">
            <span>库文件（.sdf / .smi / .csv）</span>
            <Tooltip text={PARAM_HELP.batchLibrary} />
          </div>
          <input
            type="file"
            accept=".sdf,.smi,.csv"
            onChange={(e) => setFile(e.target.files?.[0] ?? null)}
          />
          <div className="row">
            <input
              type="text"
              placeholder="库名称（可选）"
              value={libName}
              onChange={(e) => setLibName(e.target.value)}
            />
            <button onClick={upload} disabled={!file || busy}>
              {busy ? '解析中…' : '上传并解析'}
            </button>
          </div>
          {preparedLib && (
            <div className="summary">
              <p>
                成功 {preparedLib.ligand_ids.length} 条，失败 {preparedLib.failed.length} 条。
              </p>
              {preparedLib.failed.length > 0 && (
                <details>
                  <summary>失败明细（{preparedLib.failed.length}）</summary>
                  <table className="data-table">
                    <thead>
                      <tr>
                        <th>#</th>
                        <th>SMILES</th>
                        <th>原因</th>
                      </tr>
                    </thead>
                    <tbody>
                      {preparedLib.failed.map((f) => (
                        <tr key={f.index}>
                          <td>{f.index}</td>
                          <td className="mono wrap">{f.smiles || '—'}</td>
                          <td>{f.reason}</td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </details>
              )}
            </div>
          )}
        </>
      ) : !library ? (
        <p className="muted">加载中…</p>
      ) : library.length === 0 ? (
        <p className="muted">暂无已准备的配体。</p>
      ) : (
        <>
          <p className="muted">已选 {batchLigandIds.length} 个配体</p>
          <ul className="pick-list">
            {library.map((l) => (
              <li key={l.id}>
                <label className="pick-item check-item">
                  <input
                    type="checkbox"
                    checked={batchLigandIds.includes(l.id)}
                    onChange={(e) => toggle(l.id, e.target.checked)}
                  />
                  <span className="mono">{(l.name ?? l.id).slice(0, 24)}</span>
                  <span className="muted">可旋转键 {l.n_rotatable_bonds}</span>
                </label>
              </li>
            ))}
          </ul>
        </>
      )}
      {error && <p className="error-text">{error}</p>}
    </div>
  );
}
