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
}
