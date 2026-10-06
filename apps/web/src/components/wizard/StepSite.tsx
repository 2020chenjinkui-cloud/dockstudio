import { useEffect, useState } from 'react';
import { fetchFile } from '../../api';
import { HINT_POCKET_SUGGESTION, PARAM_HELP } from '../../helpText';
import type { DockBox, PreparedReceptor } from '../../types';
import MolViewer from '../../MolViewer';
import Tooltip from '../Tooltip';

interface StepSiteProps {
  receptor: PreparedReceptor;
  box: DockBox;
  onBoxChange: (box: DockBox) => void;
}

const AXES = ['x', 'y', 'z'] as const;
const AUTO_SIZE = 22.5;

export default function StepSite({ receptor, box, onBoxChange }: StepSiteProps) {
  const [receptorPdb, setReceptorPdb] = useState<string | null>(null);
  const [loadError, setLoadError] = useState<string | null>(null);

  useEffect(() => {
    setReceptorPdb(null);
    setLoadError(null);
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
            <Tooltip text={PARAM_HELP.boxCenter} />
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
            <Tooltip text={PARAM_HELP.boxSize} />
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
              onBoxChange({ center: firstCocrystal.center, size: [AUTO_SIZE, AUTO_SIZE, AUTO_SIZE] });
            }}
          >
            按共晶配体定盒（{AUTO_SIZE}³）
          </button>
          <button disabled title={HINT_POCKET_SUGGESTION}>
            按口袋建议
          </button>
        </div>
      </div>
    </div>
  );
}
