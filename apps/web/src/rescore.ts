import type { AiRescoreResult, AiRescoreTask } from './types';

export const RESCORE_TOP_DEFAULT = 3;

export function latestTaskWithResults(tasks: AiRescoreTask[]): AiRescoreTask | null {
  const withResults = tasks.filter((task) => task.results.length > 0);
  if (withResults.length === 0) return null;
  return withResults.reduce((a, b) =>
    (b.finished_at ?? b.created_at) > (a.finished_at ?? a.created_at) ? b : a,
  );
}

export function mergeRescoreResults(
  tasks: AiRescoreTask[],
): Record<number, AiRescoreResult> {
  const task = latestTaskWithResults(tasks);
  if (!task) return {};
  return Object.fromEntries(task.results.map((r) => [r.pose_index, r]));
}

export function defaultPoseSelection(nPoses: number, top: number = RESCORE_TOP_DEFAULT): number[] {
  const count = Math.max(0, Math.min(top, nPoses));
  return Array.from({ length: count }, (_, i) => i + 1);
}

export interface RescoreFormValidation {
  ok: boolean;
  error: string | null;
}

export function validateRescoreSelection(indices: number[]): RescoreFormValidation {
  if (indices.length === 0) {
    return { ok: false, error: '至少需要勾选一个姿态' };
  }
  return { ok: true, error: null };
}
