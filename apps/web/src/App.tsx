import { useCallback, useEffect, useState } from 'react';
import {
  createProject,
  deleteProject,
  getEngines,
  getHealth,
  listProjects,
} from './api';
import OnboardingOverlay, { ONBOARDING_KEY } from './components/OnboardingOverlay';
import ProjectBar from './components/ProjectBar';
import Wizard from './components/wizard/Wizard';
import DebugPage from './DebugPage';
import { useRoute } from './router';
import type { EngineManifest, Project } from './types';

export default function App() {
  const route = useRoute();
  if (route.startsWith('/debug')) {
    return <DebugPage />;
  }
  return <WizardPage />;
}

function WizardPage() {
  const [projects, setProjects] = useState<Project[] | null>(null);
  const [currentId, setCurrentId] = useState<string | null>(null);
  const [engines, setEngines] = useState<EngineManifest[]>([]);
  const [backendOk, setBackendOk] = useState(false);
  const [resumeJobId, setResumeJobId] = useState<string | null>(null);
  const [resumeBatchId, setResumeBatchId] = useState<string | null>(null);
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
    getHealth()
      .then(() => setBackendOk(true))
      .catch(() => setBackendOk(false));
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

  return (
    <div className="app">
      <ProjectBar
        projects={projects ?? []}
        current={current}
        backendOk={backendOk}
        onSelect={(id) => {
          setCurrentId(id);
          setResumeJobId(null);
          setResumeBatchId(null);
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
            {createError && <p className="error-text">{createError}</p>}
          </section>
        ) : (
          <Wizard
            key={current.id}
            projectId={current.id}
            engines={engines}
            resumeJobId={resumeJobId}
            resumeBatchId={resumeBatchId}
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
      </main>
      {showOnboarding && <OnboardingOverlay onDone={closeOnboarding} />}
    </div>
  );
}
