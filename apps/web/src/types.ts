export interface EngineManifest {
  name: string;
  version: string;
  scoring_functions: string[];
  supports_flex_residues: boolean;
  supports_batch_ligands: boolean;
  requires_gpu: boolean;
  platforms: string[];
  homepage: string | null;
  available: boolean;
  unavailable_reason: string | null;
}

export interface CocrystalLigand {
  resname: string;
  chain: string;
  n_atoms: number;
  center: [number, number, number];
}

export interface PreparedReceptor {
  id: string;
  pdbqt_path: string;
  original_path: string;
  chains: string[];
  cocrystal_ligands: CocrystalLigand[];
  warnings: string[];
  parent_id: string | null;
  flex_residues: string[];
  flex_pdbqt_path: string | null;
  receptor_json_path: string | null;
}

export interface PreparedLigand {
  id: string;
  source: 'file' | 'smiles';
  pdbqt_path: string;
  n_rotatable_bonds: number;
  smiles: string | null;
  name: string | null;
}

export interface DockBox {
  center: [number, number, number];
  size: [number, number, number];
}

export interface DockParams {
  exhaustiveness: number;
  num_modes: number;
  energy_range: number;
  cpu: number;
  seed: number | null;
  scoring: string;
}

export interface DockPose {
  index: number;
  affinity: number;
  rmsd_lb: number;
  rmsd_ub: number;
  file: string;
}

export interface DockResult {
  engine: string;
  engine_version: string;
  params: DockParams;
  box: DockBox;
  poses: DockPose[];
  log: string;
  run_dir: string;
}

export interface DockRequest {
  engine: string;
  receptor_id: string;
  ligand_id: string;
  box: DockBox;
  params: DockParams;
}

export interface HealthStatus {
  status: string;
  version: string;
}

// ---- API v0.2 ----

export interface Project {
  id: string;
  name: string;
  path: string;
  created_at: string;
}

export type JobStatus =
  | 'queued'
  | 'running'
  | 'completed'
  | 'failed'
  | 'canceled'
  | 'interrupted';

export interface ProjectJobSummary {
  id: string;
  status: JobStatus;
  ligand_id: string;
  engine: string;
  top_affinity: number | null;
  created_at: string;
  finished_at: string | null;
  batch_id?: string | null;
}

export interface ProjectDetail extends Project {
  jobs: ProjectJobSummary[];
}

export interface JobProgress {
  stage: string;
  percent: number;
  message: string;
}

export interface Job {
  id: string;
  project_id: string;
  batch_id: string | null;
  status: JobStatus;
  progress: JobProgress;
  request: DockRequest;
  result: DockResult | null;
  error: string | null;
  created_at: string;
  started_at: string | null;
  finished_at: string | null;
}

export interface CreateJobResponse {
  job_id: string;
}

export type JobEventType = 'status' | 'progress' | 'finished';

export interface JobEvent {
  type: JobEventType;
  data: unknown;
}

export interface Interaction {
  type: string;
  residue: string;
  distance: number;
  ligand_atom_ids: number[];
  detail: string;
}

export interface PoseAnalysis {
  pose_index: number;
  source: string;
  interactions: Interaction[];
  warnings: string[];
  pose_check: PoseCheck | null;
}

// ---- API v0.4 ----

export type Severity = 'info' | 'warning' | 'error';

export interface PoseIssue {
  rule: string;
  severity: Severity;
  message: string;
}

export interface PoseCheck {
  pose_index: number;
  ok: boolean;
  source: string;
  issues: PoseIssue[];
}

export const REPORT_SECTIONS = [
  'snapshot',
  'scores',
  'interactions',
  'diagram',
  'viewer3d',
  'disclaimer',
] as const;

export type ReportSection = (typeof REPORT_SECTIONS)[number];

export interface ComparisonRow {
  job_id: string;
  ligand_id: string;
  ligand_name: string;
  smiles: string;
  engine: string;
  scoring: string;
  exhaustiveness: number;
  top_affinity: number | null;
  n_poses: number;
  status: JobStatus;
}

export interface ComparisonResult {
  rows: ComparisonRow[];
}

export interface Stats {
  projects: number;
  jobs: number;
  batches: number;
  completed_jobs: number;
  avg_dock_seconds: number | null;
  by_engine: Record<string, number>;
}

// ---- API v0.5 ----

export interface AiAvailability {
  available: boolean;
  reason: string;
  mode: 'venv' | 'cmd' | null;
}

export interface Capabilities {
  ai_rescore: AiAvailability;
}

export interface AiRescoreRequest {
  pose_indices?: number[];
  use_msa_server?: boolean;
}

export interface CreateRescoreResponse {
  task_id: string;
}

export interface AiRescoreResult {
  pose_index: number;
  log_ic50: number;
  binder_probability: number;
  rank_by_score: number;
}

export interface AiRescoreTask {
  id: string;
  job_id: string;
  status: JobStatus;
  progress: JobProgress;
  pose_indices: number[];
  results: AiRescoreResult[];
  error: string | null;
  created_at: string;
  finished_at: string | null;
  msa_server_used: boolean;
}

// ---- API v0.3 ----

export interface Pocket {
  id: string;
  center: [number, number, number];
  size: [number, number, number];
  score: number;
  n_points: number;
}

export interface PocketResponse {
  receptor_id: string;
  method: string;
  pockets: Pocket[];
}

export interface FlexReceptorRequest {
  flex_residues: string[];
}

export interface LibraryFailure {
  index: number;
  smiles: string;
  reason: string;
}

export interface PreparedLibrary {
  id: string;
  name: string;
  path: string;
  ligand_ids: string[];
  failed: LibraryFailure[];
}

export interface BatchCreateRequest {
  name?: string;
  engine: string;
  receptor_id: string;
  ligand_ids: string[];
  box: DockBox;
  params: DockParams;
}

export interface CreateBatchResponse {
  batch_id: string;
}

export interface BatchSummary {
  id: string;
  project_id: string;
  name: string;
  total: number;
  done: number;
  failed: number;
  running: number;
  queued: number;
  canceled: number;
  receptor_id: string;
  engine: string;
  created_at: string;
  is_terminal: boolean;
}

export type BatchRowStatus = JobStatus;

export interface BatchRow {
  ligand_id: string;
  ligand_name: string;
  smiles: string;
  job_id: string;
  status: BatchRowStatus;
  top_affinity: number | null;
  error: string | null;
}

export interface BatchResults {
  total: number;
  rows: BatchRow[];
}

export type BatchEventType = 'batch_progress' | 'finished';

export interface BatchEvent {
  type: BatchEventType;
  data: unknown;
}
