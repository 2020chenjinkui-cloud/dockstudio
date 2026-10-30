import { useCallback, useEffect, useState } from 'react';
import { getBatch, getJob } from '../../api';
import { defaultDockCpu } from '../../defaults';
import { useI18n } from '../../i18n';
import type {
  BatchSummary,
  DemoExample,
  DockBox,
  DockParams,
  EngineManifest,
  Job,
  PreparedLigand,
  PreparedReceptor,
} from '../../types';
import BatchPanel from './BatchPanel';
import StepBatch from './StepBatch';
import StepLigand from './StepLigand';
import StepReceptor from './StepReceptor';
import StepResults from './StepResults';
import StepRun from './StepRun';
import StepSite from './StepSite';

export type WizardMode = 'single' | 'batch';

export interface WizardData {
  mode: WizardMode;
  engine: string;
  receptor: PreparedReceptor | null;
  ligand: PreparedLigand | null;
  batchLigandIds: string[];
  box: DockBox;
  params: DockParams;
  flexText: string;
  jobId: string | null;
  batchId: string | null;
}

const DEFAULT_DATA: WizardData = {
  mode: 'single',
  engine: 'vina',
  receptor: null,
  ligand: null,
  batchLigandIds: [],
  jobId: null,
  batchId: null,
  flexText: '',
  box: { center: [0, 0, 0], size: [25, 25, 25] },
  params: {
    exhaustiveness: 8,
    num_modes: 9,
    energy_range: 3.0,
    cpu: defaultDockCpu(),
    seed: null,
    scoring: 'vina',
  },
};

const STEP_KEYS = ['step.receptor', 'step.ligand', 'step.site', 'step.run', 'step.results'];

interface WizardProps {
  projectId: string;
  engines: EngineManifest[];
  resumeJobId: string | null;
  resumeBatchId: string | null;
  demo?: DemoExample | null;
  onResumeHandled: () => void;
  onOpenJob: (jobId: string) => void;
}

export default function Wizard({
  projectId,
  engines,
  resumeJobId,
  resumeBatchId,
  demo = null,
  onResumeHandled,
  onOpenJob,
}: WizardProps) {
  const { t } = useI18n();
  const [step, setStep] = useState(1);
  const [data, setData] = useState<WizardData>(DEFAULT_DATA);
  const [job, setJob] = useState<Job | null>(null);
  const [batchSummary, setBatchSummary] = useState<BatchSummary | null>(null);
  const [advanced, setAdvanced] = useState(false);
  const [resumeError, setResumeError] = useState<string | null>(null);
  const [demoNote, setDemoNote] = useState<string | null>(null);

  useEffect(() => {
    if (!demo) return;
    setData((d) => ({
      ...d,
      receptor: demo.receptor,
      ligand: demo.ligand,
      box: demo.box,
    }));
    setDemoNote(demo.note);
    setStep(3);
  }, [demo]);

  const patch = useCallback((p: Partial<WizardData>) => {
    setData((d) => ({ ...d, ...p }));
  }, []);

  const onJobUpdate = useCallback((j: Job) => setJob(j), []);

  const resetRun = useCallback(() => {
    patch({ jobId: null, batchId: null });
    setJob(null);
    setBatchSummary(null);
    setStep(4);
  }, [patch]);

  useEffect(() => {
    if (!resumeJobId) return;
    setResumeError(null);
    getJob(resumeJobId)
      .then((j) => {
        setJob(j);
        setData((d) => ({ ...d, mode: 'single', jobId: j.id }));
        setStep(j.status === 'completed' ? 5 : 4);
      })
      .catch((e: unknown) => setResumeError(e instanceof Error ? e.message : String(e)))
      .finally(onResumeHandled);
  }, [resumeJobId, onResumeHandled]);

  useEffect(() => {
    if (!resumeBatchId) return;
    setResumeError(null);
    getBatch(resumeBatchId)
      .then((s) => {
        setBatchSummary(s);
        setData((d) => ({ ...d, mode: 'batch', batchId: s.id }));
        setStep(s.is_terminal ? 5 : 4);
      })
      .catch((e: unknown) => setResumeError(e instanceof Error ? e.message : String(e)))
      .finally(onResumeHandled);
  }, [resumeBatchId, onResumeHandled]);

  const canEnter = (target: number): boolean => {
    if (target <= step) return true;
    switch (target) {
      case 2:
        return data.receptor !== null;
      case 3:
        return data.mode === 'single' ? data.ligand !== null : data.batchLigandIds.length > 0;
      case 4:
        return (
          data.receptor !== null &&
          (data.mode === 'single' ? data.ligand !== null : data.batchLigandIds.length > 0)
        );
      case 5:
        return data.mode === 'single'
          ? job?.status === 'completed'
          : batchSummary?.is_terminal === true;
      default:
        return false;
    }
  };

  const switchMode = (mode: WizardMode) => {
    if (mode === data.mode) return;
    setData((d) => ({ ...d, mode, jobId: null, batchId: null }));
    setJob(null);
    setBatchSummary(null);
    setStep((s) => (s > 3 ? 2 : s));
  };

  return (
    <div className="wizard">
      <ol className="step-bar">
        {STEP_KEYS.map((key, i) => {
          const n = i + 1;
          const reachable = canEnter(n);
          const text = n === 2 && data.mode === 'batch' ? t('step.ligandLibrary') : t(key);
          return (
            <li key={key}>
              <button
                className={`step ${n === step ? 'current' : ''} ${n < step ? 'done' : ''}`}
                disabled={!reachable}
                onClick={() => setStep(n)}
              >
                <span className="step-num">{n}</span>
                {text}
              </button>
            </li>
          );
        })}
        <li className="mode-slot">
          <button
            className={data.mode === 'single' ? 'chip active' : 'chip'}
            onClick={() => switchMode('single')}
          >
            {t('mode.single')}
          </button>
          <button
            className={data.mode === 'batch' ? 'chip active' : 'chip'}
            onClick={() => switchMode('batch')}
          >
            {t('mode.batch')}
          </button>
        </li>
        <li className="adv-slot">
          <button className={advanced ? 'chip active' : 'chip'} onClick={() => setAdvanced((v) => !v)}>
            {advanced ? t('adv.on') : t('adv.toggle')}
          </button>
        </li>
      </ol>
      {resumeError && <div className="banner">任务加载失败：{resumeError}</div>}
      {demoNote && (
        <div className="demo-note">
          <span>{demoNote}</span>
          <button className="link-btn" onClick={() => setDemoNote(null)}>
            {t('bar.close')}
          </button>
        </div>
      )}
      <section className="card wizard-card">
        {step === 1 && (
          <StepReceptor
            receptor={data.receptor}
            onPrepared={(receptor) => {
              patch({ receptor });
              setStep(2);
            }}
            onAutoBox={(center) => patch({ box: { ...data.box, center } })}
          />
        )}
        {step === 2 && (
          <StepLigand
            mode={data.mode}
            ligand={data.ligand}
            onPrepared={(ligand) => {
              patch({ ligand });
              setStep(3);
            }}
            batchLigandIds={data.batchLigandIds}
            onBatchChange={(batchLigandIds) => patch({ batchLigandIds })}
          />
        )}
        {step === 3 && data.receptor && (
          <StepSite receptor={data.receptor} box={data.box} onBoxChange={(box) => patch({ box })} />
        )}
        {step === 4 &&
          (data.mode === 'single' ? (
            <StepRun
              projectId={projectId}
              data={data}
              advanced={advanced}
              engines={engines}
              onParamsChange={(params) => patch({ params })}
              onEngineChange={(engine) => patch({ engine })}
              onFlexText={(flexText) => patch({ flexText })}
              onJobId={(jobId) => patch({ jobId })}
              onJobUpdate={onJobUpdate}
              onCompleted={() => setStep(5)}
              onReset={resetRun}
            />
          ) : (
            <StepBatch
              projectId={projectId}
              data={data}
              engines={engines}
              advanced={advanced}
              onBatchId={(batchId) => patch({ batchId })}
              onParamsChange={(params) => patch({ params })}
              onEngineChange={(engine) => patch({ engine })}
              onOpenJob={onOpenJob}
              onTerminal={(s) => {
                setBatchSummary(s);
                setStep(5);
              }}
            />
          ))}
        {step === 5 &&
          (data.mode === 'single' ? (
            data.jobId && <StepResults jobId={data.jobId} job={job} onJobUpdate={onJobUpdate} />
          ) : data.batchId ? (
            <div className="step-body">
              <div className="progress-head">
                <span
                  className={
                    batchSummary?.is_terminal ? 'badge status-completed' : 'badge status-running'
                  }
                >
                  {batchSummary?.is_terminal ? '已完成' : '进行中'}
                </span>
                {batchSummary && (
                  <span>
                    完成 {batchSummary.done}/{batchSummary.total}，失败 {batchSummary.failed}
                  </span>
                )}
              </div>
              <BatchPanel batchId={data.batchId} onOpenJob={onOpenJob} />
            </div>
          ) : null)}
      </section>
      <div className="wizard-nav">
        {step > 1 && <button onClick={() => setStep(step - 1)}>{t('nav.prev')}</button>}
        {step < 4 && (
          <button
            className="primary"
            disabled={!canEnter(step + 1)}
            onClick={() => setStep(step + 1)}
          >
            {t('nav.next')}
          </button>
        )}
        {step === 5 && (
          <button className="primary" onClick={resetRun}>
            {t('nav.newRun')}
          </button>
        )}
      </div>
    </div>
  );
}
