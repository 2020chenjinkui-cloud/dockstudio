import { useCallback, useEffect, useState } from 'react';
import { getJob } from '../../api';
import { HINT_ADVANCED_MODE } from '../../helpText';
import type {
  DockBox,
  DockParams,
  EngineManifest,
  Job,
  PreparedLigand,
  PreparedReceptor,
} from '../../types';
import StepLigand from './StepLigand';
import StepReceptor from './StepReceptor';
import StepResults from './StepResults';
import StepRun from './StepRun';
import StepSite from './StepSite';

export interface WizardData {
  receptor: PreparedReceptor | null;
  ligand: PreparedLigand | null;
  box: DockBox;
  params: DockParams;
  jobId: string | null;
}

const DEFAULT_DATA: WizardData = {
  receptor: null,
  ligand: null,
  jobId: null,
  box: { center: [0, 0, 0], size: [25, 25, 25] },
  params: {
    exhaustiveness: 8,
    num_modes: 9,
    energy_range: 3.0,
    cpu: 1,
    seed: null,
    scoring: 'vina',
  },
};

const STEPS = ['受体', '配体', '对接位点', '运行', '结果'];

interface WizardProps {
  projectId: string;
  engines: EngineManifest[];
  resumeJobId: string | null;
  onResumeHandled: () => void;
}

export default function Wizard({ projectId, engines, resumeJobId, onResumeHandled }: WizardProps) {
  const [step, setStep] = useState(1);
  const [data, setData] = useState<WizardData>(DEFAULT_DATA);
  const [job, setJob] = useState<Job | null>(null);
  const [resumeError, setResumeError] = useState<string | null>(null);

  const patch = useCallback((p: Partial<WizardData>) => {
    setData((d) => ({ ...d, ...p }));
  }, []);

  const onJobUpdate = useCallback((j: Job) => setJob(j), []);

  useEffect(() => {
    if (!resumeJobId) return;
    setResumeError(null);
    getJob(resumeJobId)
      .then((j) => {
        setJob(j);
        setData((d) => ({ ...d, jobId: j.id }));
        setStep(j.status === 'completed' ? 5 : 4);
      })
      .catch((e: unknown) => setResumeError(e instanceof Error ? e.message : String(e)))
      .finally(onResumeHandled);
  }, [resumeJobId, onResumeHandled]);

  const canEnter = (target: number): boolean => {
    if (target <= step) return true;
    switch (target) {
      case 2:
        return data.receptor !== null;
      case 3:
        return data.ligand !== null;
      case 4:
        return data.receptor !== null && data.ligand !== null;
      case 5:
        return job?.status === 'completed';
      default:
        return false;
    }
  };

  return (
    <div className="wizard">
      <ol className="step-bar">
        {STEPS.map((label, i) => {
          const n = i + 1;
          const reachable = canEnter(n);
          return (
            <li key={label}>
              <button
                className={`step ${n === step ? 'current' : ''} ${n < step ? 'done' : ''}`}
                disabled={!reachable}
                onClick={() => setStep(n)}
              >
                <span className="step-num">{n}</span>
                {label}
              </button>
            </li>
          );
        })}
        <li className="adv-slot">
          <span className="adv-link" title={HINT_ADVANCED_MODE}>
            高级模式
          </span>
        </li>
      </ol>
      {resumeError && <div className="banner">任务加载失败：{resumeError}</div>}
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
            ligand={data.ligand}
            onPrepared={(ligand) => {
              patch({ ligand });
              setStep(3);
            }}
          />
        )}
        {step === 3 && data.receptor && (
          <StepSite receptor={data.receptor} box={data.box} onBoxChange={(box) => patch({ box })} />
        )}
        {step === 4 && (
          <StepRun
            projectId={projectId}
            data={data}
            engines={engines}
            onParamsChange={(params) => patch({ params })}
            onJobId={(jobId) => patch({ jobId })}
            onJobUpdate={onJobUpdate}
            onCompleted={() => setStep(5)}
          />
        )}
        {step === 5 && data.jobId && (
          <StepResults jobId={data.jobId} job={job} onJobUpdate={onJobUpdate} />
        )}
      </section>
      <div className="wizard-nav">
        {step > 1 && (
          <button onClick={() => setStep(step - 1)}>上一步</button>
        )}
        {step < 4 && (
          <button className="primary" disabled={!canEnter(step + 1)} onClick={() => setStep(step + 1)}>
            下一步
          </button>
        )}
      </div>
    </div>
  );
}
