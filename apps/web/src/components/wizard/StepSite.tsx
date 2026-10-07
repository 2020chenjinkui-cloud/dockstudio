import { useEffect, useState } from 'react';
import { fetchFile, predictPockets } from '../../api';
import { useI18n } from '../../i18n';
import MolViewer from '../../MolViewer';
import type { DockBox, Pocket, PreparedReceptor } from '../../types';
import Tooltip from '../Tooltip';

interface StepSiteProps {
  receptor: PreparedReceptor;
  box: DockBox;
  onBoxChange: (box: DockBox) => void;
}

const AXES = ['x', 'y', 'z'] as const;
const AUTO_SIZE = 22.5;

export default function StepSite({ receptor, box, onBoxChange }: StepSiteProps) {
  const { t } = useI18n();
  const [receptorPdb, setReceptorPdb] = useState<string | null>(null);
  const [loadError, setLoadError] = useState<string | null>(null);
  const [pockets, setPockets] = useState<Pocket[] | null>(null);
  const [pocketBusy, setPocketBusy] = useState(false);
  const [pocketError, setPocketError] = useState<string | null>(null);
  const [selectedPocket, setSelectedPocket] = useState<string | null>(null);

  useEffect(() => {
    setReceptorPdb(null);
    setLoadError(null);
    setPockets(null);
    setSelectedPocket(null);
    setPocketError(null);
    let cancelled = false;
    fetchFile(receptor.original_path)
      .then((text) => {
        if (!cancelled) setReceptorPdb(text);
      })
      .catch((e: unknown) => {
        if (!cancelled) setLoadError(e instanceof Error ? e.message : String(e));
      });
    return () => {
      cancelled = true;
    };
  }, [receptor]);

  const setVec = (key: 'center' | 'size', i: number, v: number) => {
    const next: [number, number, number] = [...box[key]];
    next[i] = Number.isNaN(v) ? 0 : v;
    onBoxChange({ ...box, [key]: next });
  };

  const firstCocrystal = receptor.cocrystal_ligands[0] ?? null;

  const predict = async () => {
    setPocketBusy(true);
    setPocketError(null);
    try {
      const res = await predictPockets(receptor.id, 5);
      const sorted = [...res.pockets].sort((a, b) => b.score - a.score);
      setPockets(sorted);
      setSelectedPocket(null);
      if (sorted.length === 0) setPocketError('未预测到可用口袋，请手动定盒。');
    } catch (e) {
      setPocketError(e instanceof Error ? e.message : String(e));
    } finally {
      setPocketBusy(false);
    }
  };

  const maxScore = pockets && pockets.length > 0 ? Math.max(...pockets.map((p) => p.score)) : 1;

  return (
    <div className="step-body">
      <div className="site-layout">
        <div>
          <MolViewer
            receptorPdb={receptorPdb}
            poses={[]}
            selectedIndex={null}
            boxCenter={box.center}
            boxSize={box.size}
          />
          {loadError && <p className="error-text">受体坐标加载失败：{loadError}</p>}
        </div>
        <div className="site-controls">
          <div className="control-row">
            <span>中心（Å）</span>
            <Tooltip text={t('help.boxCenter')} />
          </div>
          <div className="triple">
            {AXES.map((axis, i) => (
              <label key={axis}>
                {axis}
                <input
                  type="number"
                  step={0.1}
                  value={box.center[i]}
                  onChange={(e) => setVec('center', i, Number.parseFloat(e.target.value))}
                />
              </label>
            ))}
          </div>
          <div className="control-row">
            <span>尺寸（Å）</span>
            <Tooltip text={t('help.boxSize')} />
          </div>
          <div className="triple">
            {AXES.map((axis, i) => (
              <label key={axis}>
                {axis}
                <input
                  type="number"
                  step={0.5}
                  min={1}
                  value={box.size[i]}
                  onChange={(e) => setVec('size', i, Number.parseFloat(e.target.value))}
                />
              </label>
            ))}
          </div>
          <button
            disabled={!firstCocrystal}
            title={firstCocrystal ? undefined : '受体中没有共晶配体'}
            onClick={() => {
              if (!firstCocrystal) return;
              setSelectedPocket(null);
              onBoxChange({ center: firstCocrystal.center, size: [AUTO_SIZE, AUTO_SIZE, AUTO_SIZE] });
            }}
          >
            按共晶配体定盒（{AUTO_SIZE}³）
          </button>
          <button onClick={() => void predict()} disabled={pocketBusy}>
            {pocketBusy ? '预测中…' : '按口袋建议'}
          </button>
          {pocketError && <p className="error-text">{pocketError}</p>}
        </div>
      </div>
      {pockets && pockets.length > 0 && (
        <ul className="pocket-list">
          {pockets.map((p) => (
            <li key={p.id}>
              <button
                className={selectedPocket === p.id ? 'pocket-item selected' : 'pocket-item'}
                onClick={() => {
                  setSelectedPocket(p.id);
                  onBoxChange({ center: p.center, size: p.size });
                }}
              >
                <span className="mono pocket-id">{p.id}</span>
                <span className="pocket-bar">
                  <span
                    className="pocket-bar-fill"
                    style={{ width: `${Math.max(4, (p.score / maxScore) * 100)}%` }}
                  />
                </span>
                <span className="mono">{p.score.toFixed(2)}</span>
                <span className="muted">
                  中心 ({p.center.map((v) => v.toFixed(1)).join(', ')})，尺寸{' '}
                  {p.size.map((v) => v.toFixed(1)).join(' × ')}
                </span>
              </button>
            </li>
          ))}
        </ul>
      )}
    </div>
  );
}
