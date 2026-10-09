import { describe, expect, it } from 'vitest';
import {
  defaultPoseSelection,
  latestTaskWithResults,
  mergeRescoreResults,
  validateRescoreSelection,
} from './rescore';
import type { AiRescoreTask } from './types';

function task(partial: Partial<AiRescoreTask>): AiRescoreTask {
  return {
    id: 't',
    job_id: 'j',
    status: 'completed',
    progress: { stage: 'finished', percent: 100, message: '' },
    pose_indices: [1, 2],
    results: [],
    error: null,
    created_at: '2026-10-08T10:00:00Z',
    finished_at: null,
    msa_server_used: false,
    ...partial,
  };
}

describe('latestTaskWithResults', () => {
  it('无任务或无结果时返回 null', () => {
    expect(latestTaskWithResults([])).toBeNull();
    expect(latestTaskWithResults([task({ id: 'a' })])).toBeNull();
  });

  it('忽略无结果的任务；按 finished_at 取最新', () => {
    const older = task({
      id: 'older',
      finished_at: '2026-10-08T10:05:00Z',
      results: [{ pose_index: 1, log_ic50: -4.1, binder_probability: 0.5, rank_by_score: 1 }],
    });
    const newer = task({
      id: 'newer',
      finished_at: '2026-10-08T11:00:00Z',
      created_at: '2026-10-08T10:50:00Z',
      results: [{ pose_index: 2, log_ic50: -5.2, binder_probability: 0.9, rank_by_score: 1 }],
    });
    const noResults = task({ id: 'empty', finished_at: '2026-10-08T12:00:00Z' });
    expect(latestTaskWithResults([older, noResults, newer])?.id).toBe('newer');
  });

  it('finished_at 缺失时回退 created_at', () => {
    const a = task({
      id: 'a',
      created_at: '2026-10-08T09:00:00Z',
      results: [{ pose_index: 1, log_ic50: -4, binder_probability: 0.4, rank_by_score: 1 }],
    });
    const b = task({
      id: 'b',
      created_at: '2026-10-08T10:00:00Z',
      results: [{ pose_index: 1, log_ic50: -5, binder_probability: 0.6, rank_by_score: 1 }],
    });
    expect(latestTaskWithResults([a, b])?.id).toBe('b');
  });
});

describe('mergeRescoreResults', () => {
  it('按 pose_index 映射权威任务结果', () => {
    const merged = mergeRescoreResults([
      task({
        results: [
          { pose_index: 1, log_ic50: -4.111, binder_probability: 0.42, rank_by_score: 2 },
          { pose_index: 3, log_ic50: -6.534, binder_probability: 0.91, rank_by_score: 1 },
        ],
      }),
    ]);
    expect(Object.keys(merged)).toHaveLength(2);
    expect(merged[1].log_ic50).toBeCloseTo(-4.111);
    expect(merged[3].binder_probability).toBeCloseTo(0.91);
    expect(merged[2]).toBeUndefined();
  });

  it('空列表返回空映射', () => {
    expect(mergeRescoreResults([])).toEqual({});
  });
});

describe('defaultPoseSelection', () => {
  it('默认 top 3', () => {
    expect(defaultPoseSelection(9)).toEqual([1, 2, 3]);
  });

  it('姿态不足 3 个时全选', () => {
    expect(defaultPoseSelection(2)).toEqual([1, 2]);
    expect(defaultPoseSelection(0)).toEqual([]);
  });
});

describe('validateRescoreSelection', () => {
  it('空选择不通过', () => {
    const r = validateRescoreSelection([]);
    expect(r.ok).toBe(false);
    expect(r.error).not.toBeNull();
  });

  it('有选择通过', () => {
    expect(validateRescoreSelection([1, 4]).ok).toBe(true);
  });
});
