export type BackendLabel = string;

const STATUS_LABELS: Record<string, string> = {
  queued: '排队中',
  running: '运行中',
  completed: '已完成',
  failed: '失败',
  canceled: '已取消',
  interrupted: '已中断',
};

const STAGE_LABELS: Record<string, string> = {
  queued: '排队等待',
  preparing: '结构准备',
  docking: '对接采样',
  scoring: '打分排序',
  exporting: '导出姿态',
  analyzing: '相互作用分析',
  finished: '完成',
  done: '完成',
};

const INTERACTION_LABELS: Record<string, string> = {
  hbond: '氢键',
  hydrophobic: '疏水作用',
  pi_stacking: 'π堆积',
  salt_bridge: '盐桥',
};

export function statusLabel(status: string): string {
  return STATUS_LABELS[status] ?? status;
}

export function stageLabel(stage: string): string {
  return STAGE_LABELS[stage] ?? stage;
}

export function interactionTypeLabel(type: string): string {
  return INTERACTION_LABELS[type] ?? type;
}

export function isTerminalStatus(status: string): boolean {
  return status === 'completed' || status === 'failed' || status === 'canceled' || status === 'interrupted';
}

export interface ResidueRef {
  resn: string;
  resi: number;
  chain?: string;
}

export function parseResidue(label: string): ResidueRef | null {
  const m = /^([A-Za-z]{1,4})\s*(-?\d+)(?::([A-Za-z0-9]+))?$/.exec(label.trim());
  if (!m) return null;
  const ref: ResidueRef = { resn: m[1], resi: Number.parseInt(m[2], 10) };
  if (m[3]) ref.chain = m[3];
  return ref;
}
