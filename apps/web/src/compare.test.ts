import { describe, expect, it } from 'vitest';
import { bestAffinity, toCompareCsv, toggleCompareId } from './compare';
import type { ComparisonRow } from './types';

describe('toggleCompareId', () => {
  it('未选则加入', () => {
    expect(toggleCompareId(['a'], 'b')).toEqual(['a', 'b']);
  });

  it('已选则移除', () => {
    expect(toggleCompareId(['a', 'b'], 'a')).toEqual(['b']);
  });
});

const ROWS: ComparisonRow[] = [
  {
    job_id: 'j1',
    ligand_id: 'l1',
    ligand_name: 'aspirin',
    smiles: 'CC(=O)Oc1ccccc1C(=O)O',
    engine: 'vina',
    scoring: 'vina',
    exhaustiveness: 8,
    top_affinity: -7.5,
    n_poses: 9,
    status: 'completed',
  },
  {
    job_id: 'j2',
    ligand_id: 'l2',
    ligand_name: 'comma,ligand',
    smiles: 'C"C',
    engine: 'smina',
    scoring: 'vina',
    exhaustiveness: 16,
    top_affinity: -8.21,
    n_poses: 5,
    status: 'completed',
  },
  {
    job_id: 'j3',
    ligand_id: 'l3',
    ligand_name: 'failed',
    smiles: 'C',
    engine: 'vina',
    scoring: 'vina',
    exhaustiveness: 8,
    top_affinity: null,
    n_poses: 0,
    status: 'failed',
  },
];

describe('bestAffinity', () => {
  it('取最小（最优）top_affinity，忽略 null', () => {
    expect(bestAffinity(ROWS)).toBe(-8.21);
  });

  it('全 null 返回 null', () => {
    expect(bestAffinity([{ ...ROWS[2] }])).toBeNull();
  });
});

describe('toCompareCsv', () => {
  it('包含表头与全部行', () => {
    const csv = toCompareCsv(ROWS);
    const lines = csv.split('\n');
    expect(lines).toHaveLength(4);
    expect(lines[0]).toBe(
      'job_id,ligand_name,smiles,engine,scoring,exhaustiveness,top_affinity,n_poses,status',
    );
    expect(lines[1]).toContain('j1');
  });

  it('逗号与引号正确转义', () => {
    const csv = toCompareCsv(ROWS);
    expect(csv).toContain('"comma,ligand"');
    expect(csv).toContain('"C""C"');
  });

  it('null affinity 输出空字段', () => {
    const csv = toCompareCsv(ROWS);
    const line3 = csv.split('\n')[3];
    expect(line3).toContain('j3');
    expect(line3).toContain(',,0,failed');
  });
});
