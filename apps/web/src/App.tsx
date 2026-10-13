import { useCallback, useEffect, useMemo, useState } from 'react';
import {
  createProject,
  deleteProject,
  getDemoExample,
  getEngines,
  getHealth,
  listProjects,
} from './api';
import BackendBootGate from './components/BackendBootGate';
import CompareDialog from './components/CompareDialog';
import OnboardingOverlay, { ONBOARDING_KEY } from './components/OnboardingOverlay';
import ProjectBar from './components/ProjectBar';
import Wizard from './components/wizard/Wizard';
import DebugPage from './DebugPage';
import { DEMO_PROJECT_NAME } from './demo';
import { I18nContext, LANG_KEY, translate, useI18n, type Lang } from './i18n';
import { useRoute } from './router';
import type { DemoExample, EngineManifest, Project } from './types';

export default function App() {
  const route = useRoute();
  const [lang, setLangState] = useState<Lang>(() => {
    const saved = localStorage.getItem(LANG_KEY);
    return saved === 'en' ? 'en' : 'zh';
  });

  const setLang = useCallback((next: Lang) => {
    localStorage.setItem(LANG_KEY, next);
    setLangState(next);
  }, []);

  const i18n = useMemo(
    () => ({ lang, setLang, t: (key: string, fallback?: string) => translate(lang, key, fallback) }),
    [lang, setLang],
  );

  return (
    <I18nContext.Provider value={i18n}>
      {route.startsWith('/debug') ? <DebugPage /> : <WizardPage />}
    </I18nContext.Provider>
  );
}

function WizardPage() {
  const { t } = useI18n();
  const [projects, setProjects] = useState<Project[] | null>(null);
  const [currentId, setCurrentId] = useState<string | null>(null);
  const [engines, setEngines] = useState<EngineManifest[]>([]);
  const [backendOk, setBackendOk] = useState(false);
  const [resumeJobId, setResumeJobId] = useState<string | null>(null);
  const [resumeBatchId, setResumeBatchId] = useState<string | null>(null);
  const [compareIds, setCompareIds] = useState<string[]>([]);
  const [compareOpen, setCompareOpen] = useState(false);
  const [demo, setDemo] = useState<DemoExample | null>(null);
  const [demoBusy, setDemoBusy] = useState(false);
  const [newName, setNewName] = useState('');
  const [createError, setCreateError] = useState<string | null>(null);
  const [showOnboarding, setShowOnboarding] = useState(
    () => localStorage.getItem(ONBOARDING_KEY) !== '1',
  );

  const refreshProjects = useCallback(async (preferId?: string) => {
    const list = await listProjects();
    setProjects(list);
    setCurrentId((cur) => {
      if (preferId && list.some((p) => p.id === preferId)) return preferId;
      if (cur && list.some((p) => p.id === cur)) return cur;
      return list[0]?.id ?? null;
    });
  }, []);

  useEffect(() => {
    // 周期心跳：一次性探测撞上瞬态失败会把状态灯永久卡住，必须轮询恢复
    const check = () => {
      getHealth()
        .then(() => setBackendOk(true))
        .catch(() => setBackendOk(false));
    };
    check();
    const timer = window.setInterval(check, 5000);
    return () => window.clearInterval(timer);
  }, []);

  useEffect(() => {
    getEngines()
      .then(setEngines)
      .catch(() => setEngines([]));
    refreshProjects().catch(() => setProjects([]));
  }, [refreshProjects]);

  const handleCreate = useCallback(
    async (name: string) => {
      const p = await createProject(name);
      await refreshProjects(p.id);
    },
    [refreshProjects],
  );

  const handleDelete = useCallback(async () => {
    if (!currentId) return;
    await deleteProject(currentId);
    setResumeJobId(null);
    setResumeBatchId(null);
    setCompareIds([]);
    await refreshProjects();
  }, [currentId, refreshProjects]);

  const closeOnboarding = useCallback(() => {
    localStorage.setItem(ONBOARDING_KEY, '1');
    setShowOnboarding(false);
  }, []);

  const current = projects?.find((p) => p.id === currentId) ?? null;

  const submitFirstProject = async () => {
    const name = newName.trim();
    if (!name) return;
    setCreateError(null);
    try {
      await handleCreate(name);
    } catch (e) {
      setCreateError(e instanceof Error ? e.message : String(e));
    }
  };

  const loadDemoProject = async () => {
    setDemoBusy(true);
    setCreateError(null);
    try {
      const demoData = await getDemoExample();
      const existing = (projects ?? []).find((p) => p.name === DEMO_PROJECT_NAME);
      const project = existing ?? (await createProject(DEMO_PROJECT_NAME));
      setDemo(demoData);
      await refreshProjects(project.id);
    } catch (e) {
      setCreateError(e instanceof Error ? e.message : String(e));
    } finally {
      setDemoBusy(false);
    }
  };

  return (
    <div className="app">
      <ProjectBar
        projects={projects ?? []}
        current={current}
        backendOk={backendOk}
        compareIds={compareIds}
        onCompareIdsChange={setCompareIds}
        onOpenCompare={() => setCompareOpen(true)}
        onSelect={(id) => {
          setCurrentId(id);
          setResumeJobId(null);
          setResumeBatchId(null);
          setCompareIds([]);
          setDemo(null);
        }}
        onCreate={handleCreate}
        onDelete={handleDelete}
        onOpenJob={(jobId) => {
          setResumeBatchId(null);
          setResumeJobId(jobId);
        }}
        onOpenBatch={(batchId) => {
          setResumeJobId(null);
          setResumeBatchId(batchId);
        }}
      />
      <main className="main-area">
        <BackendBootGate>
          {projects === null ? (
            <p className="muted">加载项目列表…</p>
          ) : projects.length === 0 || !current ? (
            <section className="card empty-state">
              <h2>欢迎使用 DockStudio</h2>
              <p className="muted">
                从创建一个项目开始：项目用于组织你的受体、配体与对接任务。
              </p>
              <div className="row">
                <input
                  type="text"
                  placeholder="项目名称，例如：KRAS 初筛"
                  value={newName}
                  onChange={(e) => setNewName(e.target.value)}
                  onKeyDown={(e) => {
                    if (e.key === 'Enter') void submitFirstProject();
                  }}
                />
                <button className="primary" onClick={() => void submitFirstProject()}>
                  创建项目
                </button>
              </div>
              <button className="chip" onClick={() => void loadDemoProject()} disabled={demoBusy}>
                {demoBusy ? t('demo.loading') : t('demo.load')}
              </button>
              {createError && <p className="error-text">{createError}</p>}
            </section>
          ) : (
            <Wizard
              key={current.id}
              projectId={current.id}
              engines={engines}
              resumeJobId={resumeJobId}
              resumeBatchId={resumeBatchId}
              demo={demo}
              onResumeHandled={() => {
                setResumeJobId(null);
                setResumeBatchId(null);
              }}
              onOpenJob={(jobId) => {
                setResumeBatchId(null);
                setResumeJobId(jobId);
              }}
            />
          )}
        </BackendBootGate>
      </main>
      {showOnboarding && <OnboardingOverlay onDone={closeOnboarding} />}
      {compareOpen && (
        <CompareDialog jobIds={compareIds} onClose={() => setCompareOpen(false)} />
      )}
    </div>
  );
}
