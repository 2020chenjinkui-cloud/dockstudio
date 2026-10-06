import { isTerminalStatus } from './labels';
import type {
  CreateJobResponse,
  DockRequest,
  DockResult,
  EngineManifest,
  HealthStatus,
  Job,
  JobEvent,
  JobProgress,
  PoseAnalysis,
  PreparedLigand,
  PreparedReceptor,
  Project,
  ProjectDetail,
} from './types';

export class ApiError extends Error {
  status: number;

  constructor(detail: string, status: number) {
    super(detail);
    this.status = status;
  }
}

function isNetworkError(e: unknown): boolean {
  return e instanceof TypeError;
}

async function handle<T>(res: Response): Promise<T> {
  if (!res.ok) {
    let detail = `${res.status} ${res.statusText}`;
    try {
      const data = (await res.json()) as { detail?: unknown };
      if (typeof data.detail === 'string') {
        detail = data.detail;
      } else if (data.detail != null) {
        detail = JSON.stringify(data.detail);
      }
    } catch {
      // 非 JSON 错误体，保留状态行
    }
    throw new ApiError(detail, res.status);
  }
  if (res.status === 204) return undefined as T;
  return (await res.json()) as T;
}

async function req<T>(input: string, init?: RequestInit): Promise<T> {
  let res: Response;
  try {
    res = await fetch(input, init);
  } catch (e) {
    if (isNetworkError(e)) throw new ApiError('后端未启动或无法连接（http://localhost:8000）', 0);
    throw e;
  }
  return handle<T>(res);
}

function postJson<T>(url: string, body: unknown): Promise<T> {
  return req<T>(url, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(body),
  });
}

// ---- M0 调试链路 ----

export async function getHealth(): Promise<HealthStatus> {
  return req<HealthStatus>('/api/health');
}

export async function getEngines(): Promise<EngineManifest[]> {
  return req<EngineManifest[]>('/api/engines');
}

export async function prepareReceptor(file: File): Promise<PreparedReceptor> {
  const form = new FormData();
  form.append('file', file);
  return req<PreparedReceptor>('/api/prepare/receptor', { method: 'POST', body: form });
}

export interface PrepareLigandInput {
  file?: File;
  smiles?: string;
  name?: string;
}

export async function prepareLigand(input: PrepareLigandInput): Promise<PreparedLigand> {
  const form = new FormData();
  if (input.file) form.append('file', input.file);
  if (input.smiles) form.append('smiles', input.smiles);
  if (input.name) form.append('name', input.name);
  return req<PreparedLigand>('/api/prepare/ligand', { method: 'POST', body: form });
}

export async function dock(body: DockRequest): Promise<DockResult> {
  return postJson<DockResult>('/api/dock', body);
}

export async function fetchFile(path: string): Promise<string> {
  const res = await fetch(`/api/file?path=${encodeURIComponent(path)}`);
  if (!res.ok) {
    let detail = `${res.status} ${res.statusText}`;
    try {
      const data = (await res.json()) as { detail?: unknown };
      if (typeof data.detail === 'string') detail = data.detail;
    } catch {
      // 忽略
    }
    throw new ApiError(detail, res.status);
  }
  return res.text();
}

// ---- 项目 ----

export async function createProject(name: string): Promise<Project> {
  return postJson<Project>('/api/projects', { name });
}

export async function listProjects(): Promise<Project[]> {
  return req<Project[]>('/api/projects');
}

export async function getProject(pid: string): Promise<ProjectDetail> {
  return req<ProjectDetail>(`/api/projects/${pid}`);
}

export async function deleteProject(pid: string): Promise<void> {
  return req<void>(`/api/projects/${pid}`, { method: 'DELETE' });
}

// ---- 已准备列表 ----

export async function listPreparedReceptors(): Promise<PreparedReceptor[]> {
  return req<PreparedReceptor[]>('/api/prepare/receptors');
}

export async function listPreparedLigands(): Promise<PreparedLigand[]> {
  return req<PreparedLigand[]>('/api/prepare/ligands');
}

// ---- 任务 ----

export async function createJob(pid: string, body: DockRequest): Promise<CreateJobResponse> {
  return postJson<CreateJobResponse>(`/api/projects/${pid}/jobs`, body);
}

export async function getJob(jid: string): Promise<Job> {
  return req<Job>(`/api/jobs/${jid}`);
}

export async function cancelJob(jid: string): Promise<Job> {
  return req<Job>(`/api/jobs/${jid}/cancel`, { method: 'POST' });
}

export async function getAnalysis(jid: string, poseIndex: number): Promise<PoseAnalysis> {
  return req<PoseAnalysis>(`/api/jobs/${jid}/analysis?pose_index=${poseIndex}`);
}

export function reportUrl(jid: string): string {
  return `/api/jobs/${jid}/report`;
}

export interface JobSubscription {
  onProgress?: (progress: JobProgress) => void;
  onJob?: (job: Job) => void;
  onError?: (message: string) => void;
}

export function subscribeJob(jid: string, handlers: JobSubscription): () => void {
  let stopped = false;
  let pollTimer: number | null = null;
  let ws: WebSocket | null = null;

  const stop = () => {
    stopped = true;
    if (pollTimer !== null) {
      clearTimeout(pollTimer);
      pollTimer = null;
    }
    if (ws) {
      ws.onclose = null;
      ws.close();
      ws = null;
    }
  };

  const startPolling = () => {
    if (stopped || pollTimer !== null) return;
    const tick = async () => {
      if (stopped) return;
      try {
        const job = await getJob(jid);
        handlers.onJob?.(job);
        if (isTerminalStatus(job.status)) {
          stop();
          return;
        }
      } catch (e) {
        handlers.onError?.(e instanceof Error ? e.message : String(e));
      }
      if (!stopped) pollTimer = window.setTimeout(tick, 2000);
    };
    pollTimer = window.setTimeout(tick, 2000);
  };

  try {
    const proto = window.location.protocol === 'https:' ? 'wss' : 'ws';
    ws = new WebSocket(`${proto}://${window.location.host}/ws/jobs/${jid}`);
    ws.onmessage = (ev: MessageEvent<string>) => {
      if (stopped) return;
      let evt: JobEvent;
      try {
        evt = JSON.parse(ev.data) as JobEvent;
      } catch {
        return;
      }
      if (evt.type === 'progress') {
        handlers.onProgress?.(evt.data as JobProgress);
      } else {
        handlers.onJob?.(evt.data as Job);
        if (evt.type === 'finished') stop();
      }
    };
    ws.onclose = () => {
      ws = null;
      if (!stopped) startPolling();
    };
    ws.onerror = () => {
      ws?.close();
    };
  } catch {
    startPolling();
  }

  return stop;
}
