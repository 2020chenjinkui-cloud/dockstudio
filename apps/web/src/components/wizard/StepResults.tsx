import { useEffect, useMemo, useRef, useState } from 'react';
import {
  auditZipUrl,
  diagramUrl,
  fetchCrystalLigand,
  fetchFile,
  getAnalysis,
  getCrystalRmsd,
  getJob,
  listPreparedReceptors,
  reportUrl,
} from '../../api';
import { DISCLAIMER } from '../../helpText';
import { useI18n } from '../../i18n';
import { interactionTypeLabel, parseResidue, poseCheckBadge } from '../../labels';
import MolViewer, {
  type MolViewerHandle,
  type ResidueHighlight,
  type ViewerPose,
} from '../../MolViewerLazy';
import { mergeRescoreResults } from '../../rescore';
import type {
  AiRescoreResult,
  AiRescoreTask,
  CrystalRmsd,
  Job,
  PoseAnalysis,
  PoseIssue,
  PreparedReceptor,
} from '../../types';
import AiRescoreSection from '../AiRescoreSection';
import CliCard from '../CliCard';
import CrystalRmsdCell from '../CrystalRmsdCell';
import ReportDialog from '../ReportDialog';

interface StepResultsProps {
  jobId: string;
  job: Job | null;
  onJobUpdate: (job: Job) => void;
}

type RightTab = 'interactions' | 'diagram';

export default function StepResults({ jobId, job, onJobUpdate }: StepResultsProps) {
  const { t } = useI18n();
  const [loadError, setLoadError] = useState<string | null>(null);
  const [receptorPdb, setReceptorPdb] = useState<string | null>(null);
  const [receptor, setReceptor] = useState<PreparedReceptor | null>(null);
  const [poseTexts, setPoseTexts] = useState<Record<number, string>>({});
  const [poseVisible, setPoseVisible] = useState<Record<number, boolean>>({});
  const [selectedPose, setSelectedPose] = useState<number | null>(null);
  const [analyses, setAnalyses] = useState<Record<number, PoseAnalysis>>({});
  const [analysisErrors, setAnalysisErrors] = useState<Record<number, string>>({});
  const [highlight, setHighlight] = useState<ResidueHighlight | null>(null);
  const [onlyPass, setOnlyPass] = useState(false);
  const [popoverPose, setPopoverPose] = useState<number | null>(null);
  const [rightTab, setRightTab] = useState<RightTab>('interactions');
  const [diagramError, setDiagramError] = useState(false);
  const [reportOpen, setReportOpen] = useState(false);
  const [rescoreTasks, setRescoreTasks] = useState<AiRescoreTask[]>([]);
  const [crystalSdf, setCrystalSdf] = useState<string | null>(null);
  const [crystalMap, setCrystalMap] = useState<Record<number, CrystalRmsd>>({});
  const viewerRef = useRef<MolViewerHandle>(null);

  const aiMap = useMemo(() => mergeRescoreResults(rescoreTasks), [rescoreTasks]);

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
    let cancelled = false;
    setCrystalSdf(null);
    fetchCrystalLigand(jobId)
      .then((text) => {
        if (!cancelled) setCrystalSdf(text);
      })
      .catch(() => {
        // 无共晶配体时静默
      });
    return () => {
      cancelled = true;
    };
  }, [jobId]);

  useEffect(() => {
    if (!result) return;
    let cancelled = false;
    for (const pose of result.poses) {
      getCrystalRmsd(jobId, pose.index)
        .then((r) => {
          if (cancelled) return;
          setCrystalMap((m) => ({ ...m, [pose.index]: r }));
        })
        .catch(() => {
          // 无共晶配体或 pose 缺失时该列显示 —
        });
    }
    return () => {
      cancelled = true;
    };
  }, [jobId, result]);

  useEffect(() => {
    if (!result) return;
    let cancelled = false;
    for (const pose of result.poses) {
      getAnalysis(jobId, pose.index)
        .then((a) => {
          if (cancelled) return;
          setAnalyses((m) => ({ ...m, [pose.index]: a }));
        })
        .catch((e: unknown) => {
          if (cancelled) return;
          setAnalysisErrors((m) => ({
            ...m,
            [pose.index]: e instanceof Error ? e.message : String(e),
          }));
        });
    }
    return () => {
      cancelled = true;
    };
  }, [jobId, result]);

  useEffect(() => {
    setHighlight(null);
    setPopoverPose(null);
    setDiagramError(false);
  }, [selectedPose]);

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

  const selectedAnalysis = selectedPose !== null ? (analyses[selectedPose] ?? null) : null;
  const selectedAnalysisError =
    selectedPose !== null ? (analysisErrors[selectedPose] ?? null) : null;
  const selectedAnalysisLoading =
    selectedPose !== null && !analyses[selectedPose] && !analysisErrors[selectedPose];

  const groupedInteractions = useMemo(() => {
    if (!selectedAnalysis) return [];
    const groups = new Map<string, typeof selectedAnalysis.interactions>();
    for (const it of selectedAnalysis.interactions) {
      const list = groups.get(it.type) ?? [];
      list.push(it);
      groups.set(it.type, list);
    }
    return [...groups.entries()];
  }, [selectedAnalysis]);

  const visiblePoses = useMemo(() => {
    if (!result) return [];
    if (!onlyPass) return result.poses;
    return result.poses.filter((p) => analyses[p.index]?.pose_check?.ok === true);
  }, [result, onlyPass, analyses]);

  if (loadError) return <div className="banner">{loadError}</div>;
  if (!job || !result) return <p className="muted">{t('common.loading')}</p>;

  const popoverCheck = popoverPose !== null ? analyses[popoverPose]?.pose_check : null;

  return (
    <div className="step-body results-layout">
      <div>
        <div className="cli-head">
          <h3 className="section-title">姿态（{result.engine} {result.engine_version}）</h3>
          <label className="check-row">
            <input
              type="checkbox"
              checked={onlyPass}
              onChange={(e) => setOnlyPass(e.target.checked)}
            />
            {t('qc.onlyOk')}
          </label>
        </div>
        {result.poses.length === 0 ? (
          <p className="muted">没有输出姿态。</p>
        ) : (
          <table className="data-table">
            <thead>
              <tr>
                <th></th>
                <th>#</th>
                <th>affinity</th>
                <th title={t('ai.note.logIC50')}>{t('ai.th.logIC50')}</th>
                <th title={t('ai.note.binder')}>{t('ai.th.binder')}</th>
                <th>rmsd 区间</th>
                <th title={t('crystal.note')}>{t('crystal.th')}</th>
                <th>质检</th>
              </tr>
            </thead>
            <tbody>
              {visiblePoses.map((p) => (
                <PoseRow
                  key={p.index}
                  pose={p}
                  analysis={analyses[p.index]}
                  ai={aiMap[p.index]}
                  crystal={crystalMap[p.index]}
                  selected={selectedPose === p.index}
                  visible={poseVisible[p.index] ?? true}
                  hasText={Boolean(poseTexts[p.index])}
                  popoverOpen={popoverPose === p.index}
                  onSelect={() => setSelectedPose(p.index)}
                  onToggleVisible={(v) => setPoseVisible((cur) => ({ ...cur, [p.index]: v }))}
                  onTogglePopover={() =>
                    setPopoverPose((cur) => (cur === p.index ? null : p.index))
                  }
                />
              ))}
            </tbody>
          </table>
        )}
        {onlyPass && visiblePoses.length === 0 && result.poses.length > 0 && (
          <p className="muted">当前没有质检合格的姿态。</p>
        )}
        {Object.keys(aiMap).length > 0 && (
          <p className="muted ai-mix-note" title={t('ai.note.logIC50')}>
            {t('ai.mixNotice')}
          </p>
        )}

        <div className="tabs results-tabs">
          <button
            className={rightTab === 'interactions' ? 'active' : ''}
            onClick={() => setRightTab('interactions')}
          >
            {t('tab.interactions')}
          </button>
          <button
            className={rightTab === 'diagram' ? 'active' : ''}
            onClick={() => setRightTab('diagram')}
          >
            {t('tab.diagram2d')}
          </button>
        </div>
        {rightTab === 'interactions' ? (
          <>
            {selectedAnalysisLoading && <p className="muted">分析中…</p>}
            {selectedAnalysisError && <p className="error-text">{selectedAnalysisError}</p>}
            {selectedAnalysis && groupedInteractions.length === 0 && (
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
            {selectedAnalysis && selectedAnalysis.warnings.length > 0 && (
              <ul className="warnings">
                {selectedAnalysis.warnings.map((w, i) => (
                  <li key={i}>{w}</li>
                ))}
              </ul>
            )}
          </>
        ) : (
          <div className="diagram-box">
            {selectedPose === null ? (
              <p className="muted">请选择一个姿态。</p>
            ) : diagramError ? (
              <p className="muted">{t('diagram.error')}</p>
            ) : (
              <img
                key={selectedPose}
                src={diagramUrl(jobId, selectedPose)}
                alt={`pose ${selectedPose} 2D diagram`}
                onError={() => setDiagramError(true)}
              />
            )}
          </div>
        )}

        <CliCard request={job.request} flex={receptor?.parent_id != null} />
        <AiRescoreSection jobId={jobId} poseCount={result.poses.length} onTasks={setRescoreTasks} />
        <div className="row">
          <button className="primary" onClick={() => setReportOpen(true)}>
            {t('report.download')}
          </button>
          <a className="link-btn" href={reportUrl(jobId)} download>
            默认完整报告
          </a>
          <a className="link-btn" href={auditZipUrl(jobId)} download>
            {t('audit.download')}
          </a>
        </div>
        <p className="disclaimer">{DISCLAIMER}</p>
      </div>
      <div>
        <div className="viewer-toolbar">
          <button
            onClick={() => {
              const dataUrl = viewerRef.current?.capturePng();
              if (!dataUrl) return;
              const a = document.createElement('a');
              a.href = dataUrl;
              a.download = `dock-pose-${selectedPose ?? 0}.png`;
              a.click();
            }}
          >
            {t('shot.save')}
          </button>
          {crystalSdf !== null && <span className="muted">{t('crystal.th')}</span>}
        </div>
        <MolViewer
          ref={viewerRef}
          receptorPdb={receptorPdb}
          poses={viewerPoses}
          selectedIndex={selectedPose}
          boxCenter={result.box.center}
          boxSize={result.box.size}
          highlightResidue={highlight}
          crystalSdf={crystalSdf}
        />
        <details className="log-details">
          <summary>原始 log</summary>
          <pre className="log">{result.log}</pre>
        </details>
      </div>
      {popoverPose !== null && popoverCheck && (
        <div className="popover-overlay" onClick={() => setPopoverPose(null)}>
          <div className="popover" onClick={(e) => e.stopPropagation()}>
            <div className="cli-head">
              <h4>
                {t('qc.popoverTitle')} · pose {popoverPose}
              </h4>
              <button className="link-btn" onClick={() => setPopoverPose(null)}>
                {t('bar.close')}
              </button>
            </div>
            {popoverCheck.issues.length === 0 ? (
              <p className="ok-text">{t('qc.ok')}</p>
            ) : (
              <ul className="issue-list">
                {popoverCheck.issues.map((issue: PoseIssue, i: number) => (
                  <li key={i}>
                    <span className={`badge qc-${issue.severity}`}>
                      {t(`qc.severity.${issue.severity}`, issue.severity)}
                    </span>
                    <span className="mono">{issue.rule}</span>
                    <span>{issue.message}</span>
                  </li>
                ))}
              </ul>
            )}
          </div>
        </div>
      )}
      {reportOpen && <ReportDialog jobId={jobId} onClose={() => setReportOpen(false)} />}
    </div>
  );
}

interface PoseRowProps {
  pose: { index: number; affinity: number; rmsd_lb: number; rmsd_ub: number };
  analysis: PoseAnalysis | undefined;
  ai: AiRescoreResult | undefined;
  crystal: CrystalRmsd | undefined;
  selected: boolean;
  visible: boolean;
  hasText: boolean;
  popoverOpen: boolean;
  onSelect: () => void;
  onToggleVisible: (v: boolean) => void;
  onTogglePopover: () => void;
}

function PoseRow({
  pose,
  analysis,
  ai,
  crystal,
  selected,
  visible,
  hasText,
  popoverOpen,
  onSelect,
  onToggleVisible,
  onTogglePopover,
}: PoseRowProps) {
  const { t } = useI18n();
  const badge = poseCheckBadge(analysis?.pose_check);
  return (
    <tr className={selected ? 'selected' : ''} onClick={onSelect}>
      <td onClick={(e) => e.stopPropagation()}>
        <input
          type="checkbox"
          checked={visible}
          disabled={!hasText}
          onChange={(e) => onToggleVisible(e.target.checked)}
        />
      </td>
      <td>{pose.index}</td>
      <td>{pose.affinity.toFixed(1)}</td>
      <td title={ai ? t('ai.note.logIC50') : undefined}>
        {ai ? ai.log_ic50.toFixed(3) : '—'}
      </td>
      <td title={ai ? t('ai.note.binder') : undefined}>
        {ai ? `${(ai.binder_probability * 100).toFixed(1)}%` : '—'}
      </td>
      <td>
        {pose.rmsd_lb.toFixed(1)} – {pose.rmsd_ub.toFixed(1)}
      </td>
      <td>
        <CrystalRmsdCell rmsd={crystal?.rmsd ?? null} />
      </td>
      <td onClick={(e) => e.stopPropagation()}>
        {badge === 'none' ? (
          <span className="muted">—</span>
        ) : badge === 'ok' ? (
          <button className="badge qc-ok" onClick={onTogglePopover}>
            {t('qc.ok')}
          </button>
        ) : (
          <button
            className={popoverOpen ? `badge qc-${badge} qc-open` : `badge qc-${badge}`}
            onClick={onTogglePopover}
          >
            {t('qc.issue')}
          </button>
        )}
      </td>
    </tr>
  );
}
