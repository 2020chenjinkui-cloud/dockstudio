import type { DockBox, DockParams } from '../types';

interface BoxParamsCardProps {
  box: DockBox;
  params: DockParams;
  scoringOptions: string[];
  onBoxChange: (box: DockBox) => void;
  onParamsChange: (params: DockParams) => void;
}

const AXES = ['x', 'y', 'z'] as const;

export default function BoxParamsCard({
  box,
  params,
  scoringOptions,
  onBoxChange,
  onParamsChange,
}: BoxParamsCardProps) {
  const setVec = (key: 'center' | 'size', i: number, v: number) => {
    const next: [number, number, number] = [...box[key]];
    next[i] = v;
    onBoxChange({ ...box, [key]: next });
  };

  const num = (v: string): number => {
    const n = Number.parseFloat(v);
    return Number.isNaN(n) ? 0 : n;
  };

  const scoringList = scoringOptions.includes(params.scoring)
    ? scoringOptions
    : [params.scoring, ...scoringOptions];

  return (
    <section className="card">
      <h2>3. 盒子与参数</h2>
      <div className="grid-form">
        <span className="label">中心</span>
        {AXES.map((axis, i) => (
          <label key={axis}>
            {axis}
            <input
              type="number"
              step={0.1}
              value={box.center[i]}
              onChange={(e) => setVec('center', i, num(e.target.value))}
            />
          </label>
        ))}
        <span className="label">尺寸</span>
        {AXES.map((axis, i) => (
          <label key={axis}>
            {axis}
            <input
              type="number"
              step={0.5}
              min={1}
              value={box.size[i]}
              onChange={(e) => setVec('size', i, num(e.target.value))}
            />
          </label>
        ))}
      </div>
      <div className="grid-form">
        <label>
          exhaustiveness
          <input
            type="number"
            min={1}
            value={params.exhaustiveness}
            onChange={(e) => onParamsChange({ ...params, exhaustiveness: num(e.target.value) })}
          />
        </label>
        <label>
          num_modes
          <input
            type="number"
            min={1}
            value={params.num_modes}
            onChange={(e) => onParamsChange({ ...params, num_modes: num(e.target.value) })}
          />
        </label>
        <label>
          cpu
          <input
            type="number"
            min={1}
            value={params.cpu}
            onChange={(e) => onParamsChange({ ...params, cpu: num(e.target.value) })}
          />
        </label>
        <label>
          scoring
          <select
            value={params.scoring}
            onChange={(e) => onParamsChange({ ...params, scoring: e.target.value })}
          >
            {scoringList.map((s) => (
              <option key={s} value={s}>
                {s}
              </option>
            ))}
          </select>
        </label>
      </div>
    </section>
  );
}
