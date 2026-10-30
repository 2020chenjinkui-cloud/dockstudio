import { isTerminalStatus } from './labels';
import type {
  AiRescoreRequest,
  AiRescoreTask,
  BatchCreateRequest,
  BatchEvent,
  BatchResults,
  BatchSummary,
  Capabilities,
  ComparisonResult,
  CreateBatchResponse,
  CreateRescoreResponse,
  CreateJobResponse,
  CrystalRmsd,
  DemoExample,
  DemoLibrary,
  DockRequest,
  DockResult,
  EngineManifest,
  HealthStatus,
  Job,
  JobEvent,
  JobProgress,
  PocketResponse,
  PoseAnalysis,
  PreparedLibrary,
  PreparedLigand,
  PreparedReceptor,
  Project,
  ProjectDetail,
  Stats,
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

// 桌面端（Tauri webview）里相对路径会解析到 tauri.localhost 协议，到不了后端，
// 必须显式指向 127.0.0.1:8000；浏览器/开发模式保持相对路径走 vite 代理。
export const IS_TAURI =
  typeof window !== 'undefined' &&
  (window.location.protocol === 'tauri:' ||
    window.location.hostname === 'tauri.localhost' ||
    '__TAURI_INTERNALS__' in window);
export const API_BASE = IS_TAURI ? 'http://127.0.0.1:8000' : '';
const WS_HOST = IS_TAURI ? '127.0.0.1:8000' : typeof window !== 'undefined' ? window.location.host : '';

/** Absolute API URL for links/downloads/images (Tauri webview needs the host). */
export function apiUrl(path: string): string {
  return API_BASE + path;
}

async function req<T>(input: string, init?: RequestInit): Promise<T> {
  let res: Response;
  try {
    res = await fetch(API_BASE + input, init);
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
  const res = await fetch(`${API_BASE}/api/file?path=${encodeURIComponent(path)}`);
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
  return apiUrl(`/api/jobs/${jid}/report`);
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
    ws = new WebSocket(`${proto}://${WS_HOST}/ws/jobs/${jid}`);
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

// ---- 口袋预测 ----

export async function predictPockets(receptorId: string, topN = 5): Promise<PocketResponse> {
  return postJson<PocketResponse>('/api/pockets', { receptor_id: receptorId, top_n: topN });
}

// ---- 柔性残基 ----

export async function createFlexReceptor(
  rid: string,
  flexResidues: string[],
): Promise<PreparedReceptor> {
  return postJson<PreparedReceptor>(`/api/prepare/receptor/${rid}/flex`, {
    flex_residues: flexResidues,
  });
}

// ---- 配体库 ----

export async function prepareLibrary(file: File, name?: string): Promise<PreparedLibrary> {
  const form = new FormData();
  form.append('file', file);
  if (name) form.append('name', name);
  return req<PreparedLibrary>('/api/prepare/library', { method: 'POST', body: form });
}

// ---- 批量筛选 ----

export async function createBatch(
  pid: string,
  body: BatchCreateRequest,
): Promise<CreateBatchResponse> {
  return postJson<CreateBatchResponse>(`/api/projects/${pid}/batches`, body);
}

export async function getBatch(bid: string): Promise<BatchSummary> {
  return req<BatchSummary>(`/api/batches/${bid}`);
}

export async function cancelBatch(bid: string): Promise<BatchSummary> {
  return req<BatchSummary>(`/api/batches/${bid}/cancel`, { method: 'POST' });
}

export interface BatchResultsQuery {
  sort?: string;
  order?: 'asc' | 'desc';
  limit?: number;
  offset?: number;
}

export async function getBatchResults(
  bid: string,
  query: BatchResultsQuery = {},
): Promise<BatchResults> {
  const params = new URLSearchParams();
  if (query.sort) params.set('sort', query.sort);
  if (query.order) params.set('order', query.order);
  if (query.limit !== undefined) params.set('limit', String(query.limit));
  if (query.offset !== undefined) params.set('offset', String(query.offset));
  const qs = params.toString();
  return req<BatchResults>(`/api/batches/${bid}/results${qs ? `?${qs}` : ''}`);
}

export function batchResultsCsvUrl(bid: string): string {
  return apiUrl(`/api/batches/${bid}/results.csv`);
}

export interface BatchSubscription {
  onSummary?: (summary: BatchSummary) => void;
  onError?: (message: string) => void;
}

export function subscribeBatch(bid: string, handlers: BatchSubscription): () => void {
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
        const summary = await getBatch(bid);
        handlers.onSummary?.(summary);
        if (summary.is_terminal) {
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
    ws = new WebSocket(`${proto}://${WS_HOST}/ws/batches/${bid}`);
    ws.onmessage = (ev: MessageEvent<string>) => {
      if (stopped) return;
      let evt: BatchEvent;
      try {
        evt = JSON.parse(ev.data) as BatchEvent;
      } catch {
        return;
      }
      const summary = evt.data as BatchSummary;
      handlers.onSummary?.(summary);
      if (evt.type === 'finished') stop();
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

// ---- 姿态质检 / 2D 图 ----

export function diagramUrl(jid: string, poseIndex: number): string {
  return apiUrl(`/api/jobs/${jid}/diagram?pose_index=${poseIndex}`);
}

// ---- 报告 v2 ----

export interface ReportOptions {
  sections?: string[];
  title?: string;
}

export function buildReportQuery(options: ReportOptions): string {
  const params = new URLSearchParams();
  if (options.sections && options.sections.length > 0) {
    params.set('sections', options.sections.join(','));
  }
  if (options.title && options.title.trim()) {
    params.set('title', options.title.trim());
  }
  return params.toString();
}

export function reportUrlV2(jid: string, options: ReportOptions = {}): string {
  const qs = buildReportQuery(options);
  return apiUrl(`/api/jobs/${jid}/report${qs ? `?${qs}` : ''}`);
}

// ---- 结果对比 ----

export async function compareJobs(jobIds: string[]): Promise<ComparisonResult> {
  return req<ComparisonResult>(`/api/compare?job_ids=${jobIds.map(encodeURIComponent).join(',')}`);
}

// ---- 统计 ----

export async function getStats(): Promise<Stats> {
  return req<Stats>('/api/stats');
}

// ---- AI 复核（Boltz-2）----

export async function getCapabilities(): Promise<Capabilities> {
  return req<Capabilities>('/api/capabilities');
}

export async function createAiRescore(
  jid: string,
  body: AiRescoreRequest,
): Promise<CreateRescoreResponse> {
  return postJson<CreateRescoreResponse>(`/api/jobs/${jid}/rescore/ai`, body);
}

export async function listRescoreTasks(jid: string): Promise<AiRescoreTask[]> {
  return req<AiRescoreTask[]>(`/api/jobs/${jid}/rescore`);
}

export async function getRescoreTask(taskId: string): Promise<AiRescoreTask> {
  return req<AiRescoreTask>(`/api/rescore/${taskId}`);
}

export async function cancelRescoreTask(taskId: string): Promise<AiRescoreTask> {
  return req<AiRescoreTask>(`/api/rescore/${taskId}/cancel`, { method: 'POST' });
}

// ---- 演示与共晶参考 ----

export async function getDemoExample(): Promise<DemoExample> {
  return req<DemoExample>('/api/demo/example');
}

export async function getDemoLibrary(): Promise<DemoLibrary> {
  return req<DemoLibrary>('/api/demo/library');
}

export async function fetchCrystalLigand(jid: string): Promise<string | null> {
  const res = await fetch(`${API_BASE}/api/jobs/${jid}/crystal_ligand`);
  if (res.status === 404) return null;
  if (!res.ok) return null;
  return res.text();
}

export async function getCrystalRmsd(jid: string, poseIndex: number): Promise<CrystalRmsd> {
  return req<CrystalRmsd>(`/api/jobs/${jid}/crystal_rmsd?pose_index=${poseIndex}`);
}

// ---- 溯源 ----

export function auditZipUrl(jid: string): string {
  return apiUrl(`/api/jobs/${jid}/audit.zip`);
}

export function projectEventsUrl(pid: string): string {
  return apiUrl(`/api/projects/${pid}/events.jsonl`);
}
