import { useEffect, useMemo, useState } from 'react';
import {
  fetchFile,
  getAnalysis,
  getJob,
  listPreparedReceptors,
  reportUrl,
} from '../../api';
import { DISCLAIMER } from '../../helpText';
import { interactionTypeLabel, parseResidue } from '../../labels';
import MolViewer, { type ResidueHighlight, type ViewerPose } from '../../MolViewer';
import type { Job, PoseAnalysis, PreparedReceptor } from '../../types';
import CliCard from '../CliCard';

interface StepResultsProps {
  jobId: string;
  job: Job | null;
  onJobUpdate: (job: Job) => void;
}

export default function StepResults({ jobId, job, onJobUpdate }: StepResultsProps) {
  const [loadError, setLoadError] = useState<string | null>(null);
  const [receptorPdb, setReceptorPdb] = useState<string | null>(null);
  const [receptor, setReceptor] = useState<PreparedReceptor | null>(null);
  const [poseTexts, setPoseTexts] = useState<Record<number, string>>({});
  const [poseVisible, setPoseVisible] = useState<Record<number, boolean>>({});
  const [selectedPose, setSelectedPose] = useState<number | null>(null);
  const [analysis, setAnalysis] = useState<PoseAnalysis | null>(null);
  const [analysisError, setAnalysisError] = useState<string | null>(null);
  const [analysisLoading, setAnalysisLoading] = useState(false);
  const [highlight, setHighlight] = useState<ResidueHighlight | null>(null);

  useEffect(() => {
    if (job) return;
    let cancelled = false;
    getJob(jobId)
      .then((j) => {
        if (!cancelled) onJobUpdate(j);
      })
      .catch((e: unknown) => {
        if (!cancelled) setLoadError(e instanceof Error ? e.message : String(e));
      });
    return () => {
      cancelled = true;
    };
  }, [jobId, job, onJobUpdate]);

  const result = job?.result ?? null;

  useEffect(() => {
    if (!result) return;
    let cancelled = false;
    const texts: Record<number, string> = {};
    Promise.all(
      result.poses
        .filter((p) => p.file)
        .map(async (p) => {
          try {
            texts[p.index] = await fetchFile(p.file);
          } catch {
            // 单个姿态加载失败不阻塞整体
          }
        }),
    ).then(() => {
      if (cancelled) return;
      setPoseTexts(texts);
      setPoseVisible(Object.fromEntries(result.poses.map((p) => [p.index, true])));
      setSelectedPose((cur) => cur ?? result.poses[0]?.index ?? null);
    });
    return () => {
      cancelled = true;
    };
  }, [result]);

  useEffect(() => {
    if (!job) return;
    let cancelled = false;
    listPreparedReceptors()
      .then((list) => list.find((r) => r.id === job.request.receptor_id))
      .then((rec) => {
        if (cancelled || !rec) return null;
        setReceptor(rec);
        return fetchFile(rec.original_path);
      })
      .then((text) => {
        if (!cancelled && text) setReceptorPdb(text);
      })
      .catch(() => {
        // 受体缺失时仅显示姿态
      });
    return () => {
      cancelled = true;
    };
  }, [job]);

  useEffect(() => {
    if (selectedPose === null) return;
    setAnalysis(null);
    setAnalysisError(null);
    setAnalysisLoading(true);
    setHighlight(null);
    let cancelled = false;
    getAnalysis(jobId, selectedPose)
      .then((a) => {
        if (!cancelled) setAnalysis(a);
      })
      .catch((e: unknown) => {
        if (!cancelled) setAnalysisError(e instanceof Error ? e.message : String(e));
      })
      .finally(() => {
        if (!cancelled) setAnalysisLoading(false);
      });
    return () => {
      cancelled = true;
    };
  }, [jobId, selectedPose]);

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

  const groupedInteractions = useMemo(() => {
    if (!analysis) return [];
    const groups = new Map<string, typeof analysis.interactions>();
    for (const it of analysis.interactions) {
      const list = groups.get(it.type) ?? [];
      list.push(it);
      groups.set(it.type, list);
    }
    return [...groups.entries()];
  }, [analysis]);

  if (loadError) return <div className="banner">{loadError}</div>;
  if (!job || !result) return <p className="muted">结果加载中…</p>;

  return (
    <div className="step-body results-layout">
      <div>
        <h3 className="section-title">姿态（{result.engine} {result.engine_version}）</h3>
        {result.poses.length === 0 ? (
          <p className="muted">没有输出姿态。</p>
        ) : (
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
        <h3 className="section-title">相互作用（pose {selectedPose ?? '—'}）</h3>
        {analysisLoading && <p className="muted">分析中…</p>}
        {analysisError && <p className="error-text">{analysisError}</p>}
        {analysis && groupedInteractions.length === 0 && (
          <p className="muted">该姿态未检出显著相互作用。</p>
        )}
        {groupedInteractions.map(([type, list]) => (
          <div key={type} className="interaction-group">
            <h4>
              {interactionTypeLabel(type)} <span className="muted">×{list.length}</span>
            </h4>
            <ul>
              {list.map((it, i) => (
                <li key={i}>
                  <button
                    className="interaction-item"
                    onClick={() => setHighlight(parseResidue(it.residue))}
                    title="在 3D 视图中高亮该残基"
                  >
                    <span className="mono">{it.residue}</span>
                    <span>{it.distance.toFixed(2)} Å</span>
                    {it.detail && <span className="muted">{it.detail}</span>}
                  </button>
                </li>
              ))}
            </ul>
          </div>
        ))}
        {analysis && analysis.warnings.length > 0 && (
          <ul className="warnings">
            {analysis.warnings.map((w, i) => (
              <li key={i}>{w}</li>
            ))}
          </ul>
        )}
        <CliCard request={job.request} flex={receptor?.parent_id != null} />
        <a className="btn-link" href={reportUrl(jobId)} download>
          下载报告（HTML）
        </a>
        <p className="disclaimer">{DISCLAIMER}</p>
      </div>
      <div>
        <MolViewer
          receptorPdb={receptorPdb}
          poses={viewerPoses}
          selectedIndex={selectedPose}
          boxCenter={result.box.center}
          boxSize={result.box.size}
          highlightResidue={highlight}
        />
        <details className="log-details">
          <summary>原始 log</summary>
          <pre className="log">{result.log}</pre>
        </details>
      </div>
    </div>
  );
}
