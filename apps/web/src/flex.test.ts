import { describe, expect, it } from 'vitest';
import { parseFlexInput, validateFlexResidues } from './flex';

describe('parseFlexInput', () => {
  it('按空格/逗号/分号/换行拆分', () => {
    expect(parseFlexInput('A:42 A:57')).toEqual(['A:42', 'A:57']);
    expect(parseFlexInput('A:42,A:57;B:10\nC:3')).toEqual(['A:42', 'A:57', 'B:10', 'C:3']);
    expect(parseFlexInput('A:42，A:57；B:10')).toEqual(['A:42', 'A:57', 'B:10']);
  });

  it('空输入返回空数组', () => {
    expect(parseFlexInput('  ')).toEqual([]);
  });
});

describe('validateFlexResidues', () => {
  it('接受合法条目并规范化', () => {
    const r = validateFlexResidues(['A:42', 'B:7'], ['A', 'B']);
    expect(r.errors).toEqual([]);
    expect(r.valid).toEqual(['A:42', 'B:7']);
  });

  it('格式错误', () => {
    const r = validateFlexResidues(['A-42', '42:A', 'A:'], ['A']);
    expect(r.valid).toEqual([]);
    expect(r.errors).toHaveLength(3);
    expect(r.errors[0]).toContain('格式');
  });

  it('链必须存在于受体', () => {
    const r = validateFlexResidues(['Z:1'], ['A', 'B']);
    expect(r.valid).toEqual([]);
    expect(r.errors[0]).toContain('链 Z');
    expect(r.errors[0]).toContain('A / B');
  });

  it('重复条目标记错误', () => {
    const r = validateFlexResidues(['A:42', 'A:42'], ['A']);
    expect(r.valid).toEqual(['A:42']);
    expect(r.errors[0]).toContain('重复');
  });

  it('超过 10 个报错', () => {
    const many = Array.from({ length: 11 }, (_, i) => `A:${i + 1}`);
    const r = validateFlexResidues(many, ['A']);
    expect(r.errors.some((e) => e.includes('最多'))).toBe(true);
  });

  it('受体无链信息时不校验链', () => {
    const r = validateFlexResidues(['Z:1'], []);
    expect(r.errors).toEqual([]);
    expect(r.valid).toEqual(['Z:1']);
  });
});
