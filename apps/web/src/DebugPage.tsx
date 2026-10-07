import { useEffect, useMemo, useState } from 'react';
import { dock, fetchFile, getEngines, getHealth } from './api';
import BoxParamsCard from './components/BoxParamsCard';
import CliCard from './components/CliCard';
import LigandCard from './components/LigandCard';
import ReceptorCard from './components/ReceptorCard';
import MolViewer, { type ViewerPose } from './MolViewer';
import type {
  DockBox,
  DockParams,
  DockRequest,
  DockResult,
  EngineManifest,
  PreparedLigand,
  PreparedReceptor,
} from './types';

const DEFAULT_BOX: DockBox = { center: [0, 0, 0], size: [25, 25, 25] };
const DEFAULT_PARAMS: DockParams = {
  exhaustiveness: 8,
  num_modes: 9,
  energy_range: 3.0,
  cpu: 1,
  seed: null,
  scoring: 'vina',
};

export default function DebugPage() {
  const [health, setHealth] = useState<string | null>(null);
  const [engines, setEngines] = useState<EngineManifest[]>([]);
  const [backendError, setBackendError] = useState<string | null>(null);
  const [receptor, setReceptor] = useState<PreparedReceptor | null>(null);
  const [ligand, setLigand] = useState<PreparedLigand | null>(null);
  const [box, setBox] = useState<DockBox>(DEFAULT_BOX);
  const [params, setParams] = useState<DockParams>(DEFAULT_PARAMS);
  const [docking, setDocking] = useState(false);
  const [elapsed, setElapsed] = useState<number | null>(null);
  const [result, setResult] = useState<DockResult | null>(null);
  const [lastRequest, setLastRequest] = useState<DockRequest | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [receptorPdb, setReceptorPdb] = useState<string | null>(null);
  const [poseTexts, setPoseTexts] = useState<Record<number, string>>({});
  const [poseVisible, setPoseVisible] = useState<Record<number, boolean>>({});
  const [selectedPose, setSelectedPose] = useState<number | null>(null);

  useEffect(() => {
    getHealth()
      .then((h) => setHealth(`${h.status} v${h.version}`))
      .catch(() => setHealth(null));
    getEngines()
      .then(setEngines)
      .catch((e: unknown) => setBackendError(e instanceof Error ? e.message : String(e)));
  }, []);

  useEffect(() => {
    setReceptorPdb(null);
    if (!receptor) return;
    let cancelled = false;
    fetchFile(receptor.original_path)
      .then((text) => {
        if (!cancelled) setReceptorPdb(text);
      })
      .catch(() => {
        // 受体坐标加载失败不阻塞流程
      });
    return () => {
      cancelled = true;
    };
  }, [receptor]);

  const scoringOptions = useMemo(
    () => [...new Set(engines.filter((e) => e.available).flatMap((e) => e.scoring_functions))],
    [engines],
  );

  const viewerPoses: ViewerPose[] = useMemo(() => {
    if (!result) return [];
    return result.poses
      .filter((p) => poseTexts[p.index])
      .map((p) => ({
        index: p.index,
        sdf: poseTexts[p.index],
        visible: poseVisible[p.index] ?? true,
      }));
  }, [result, poseTexts, poseVisible]);

  const runDock = async () => {
    if (!receptor || !ligand) return;
    setDocking(true);
    setError(null);
    setResult(null);
    setElapsed(null);
    setPoseTexts({});
    setSelectedPose(null);
    const t0 = performance.now();
    try {
      const engine = engines.find((e) => e.available)?.name ?? 'vina';
      const request: DockRequest = {
        engine,
        receptor_id: receptor.id,
        ligand_id: ligand.id,
        box,
        params,
      };
      const res = await dock(request);
      setLastRequest(request);
      setResult(res);
      setElapsed((performance.now() - t0) / 1000);
      const texts: Record<number, string> = {};
      const failures: string[] = [];
      await Promise.all(
        res.poses
          .filter((p) => p.file)
          .map(async (p) => {
            try {
              texts[p.index] = await fetchFile(p.file);
            } catch {
              failures.push(`pose ${p.index}`);
            }
          }),
      );
      setPoseTexts(texts);
      setPoseVisible(Object.fromEntries(res.poses.map((p) => [p.index, true])));
      setSelectedPose(res.poses[0]?.index ?? null);
      if (failures.length > 0) {
        setError(`部分姿态文件加载失败：${failures.join(', ')}`);
      }
    } catch (e) {
      setError(e instanceof Error ? e.message : String(e));
    } finally {
      setDocking(false);
    }
  };

  return (
    <div className="app">
      <header>
        <h1>DockStudio · M0 调试台</h1>
        <span className="muted">
          后端：{backendError ? '未连接' : (health ?? '—')}
        </span>
      </header>
      {backendError && <div className="banner">/api/engines 加载失败：{backendError}</div>}
      {error && <div className="banner">{error}</div>}
      <div className="layout">
        <div className="left-col">
          <ReceptorCard
            receptor={receptor}
            onPrepared={setReceptor}
            onSelectCocrystal={(center) => setBox((b) => ({ ...b, center }))}
          />
          <LigandCard ligand={ligand} onPrepared={setLigand} />
          <BoxParamsCard
            box={box}
            params={params}
            scoringOptions={scoringOptions}
            onBoxChange={setBox}
            onParamsChange={setParams}
          />
          <section className="card run-card">
            <button
              className="primary"
              onClick={runDock}
              disabled={!receptor || !ligand || docking}
            >
              {docking ? '对接中…' : '运行对接'}
            </button>
            {(!receptor || !ligand) && (
              <p className="muted">需要先完成受体（1）与配体（2）的准备。</p>
            )}
            {elapsed !== null && <span className="elapsed">用时 {elapsed.toFixed(2)} s</span>}
          </section>
        </div>
        <div className="right-col">
          <section className="card">
            <h2>姿态</h2>
            {!result && <p className="muted">运行对接后在此显示结果。</p>}
            {result && result.poses.length === 0 && <p className="muted">没有输出姿态。</p>}
            {result && result.poses.length > 0 && (
              <table className="data-table">
                <thead>
                  <tr>
                    <th></th>
                    <th>#</th>
                    <th>affinity</th>
                    <th>rmsd 区间</th>
                  </tr>
                </thead>
                <tbody>
                  {result.poses.map((p) => (
                    <tr
                      key={p.index}
                      className={selectedPose === p.index ? 'selected' : ''}
                      onClick={() => setSelectedPose(p.index)}
                    >
                      <td onClick={(e) => e.stopPropagation()}>
                        <input
                          type="checkbox"
                          checked={poseVisible[p.index] ?? true}
                          disabled={!poseTexts[p.index]}
                          onChange={(e) =>
                            setPoseVisible((v) => ({ ...v, [p.index]: e.target.checked }))
                          }
                        />
                      </td>
                      <td>{p.index}</td>
                      <td>{p.affinity.toFixed(1)}</td>
                      <td>
                        {p.rmsd_lb.toFixed(1)} – {p.rmsd_ub.toFixed(1)}
                      </td>
                    </tr>
                  ))}
                </tbody>
              </table>
            )}
          </section>
          <section className="card">
            <h2>3D 视图</h2>
            <MolViewer
              receptorPdb={receptorPdb}
              poses={viewerPoses}
              selectedIndex={selectedPose}
              boxCenter={box.center}
              boxSize={box.size}
            />
          </section>
          {result && lastRequest && (
            <section className="card">
              <CliCard request={lastRequest} flex={receptor?.parent_id != null} />
            </section>
          )}
          {result && (
            <section className="card">
              <details>
                <summary>原始 log（{result.engine} {result.engine_version}）</summary>
                <pre className="log">{result.log}</pre>
              </details>
            </section>
          )}
        </div>
      </div>
      <footer>
        {engines.length === 0 && <span className="muted">引擎列表为空（后端未连接？）</span>}
        <ul>
          {engines.map((e) => (
            <li key={e.name} title={e.available ? undefined : (e.unavailable_reason ?? undefined)}>
              <span className={e.available ? 'dot ok' : 'dot bad'} />
              {e.name} {e.version}
              {!e.available && <span className="muted">（不可用）</span>}
            </li>
          ))}
        </ul>
      </footer>
    </div>
  );
}
