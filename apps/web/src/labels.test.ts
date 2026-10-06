import { describe, expect, it } from 'vitest';
import {
  interactionTypeLabel,
  isTerminalStatus,
  parseResidue,
  stageLabel,
  statusLabel,
} from './labels';

describe('statusLabel', () => {
  it('映射已知状态为中文', () => {
    expect(statusLabel('queued')).toBe('排队中');
    expect(statusLabel('running')).toBe('运行中');
    expect(statusLabel('completed')).toBe('已完成');
    expect(statusLabel('failed')).toBe('失败');
    expect(statusLabel('canceled')).toBe('已取消');
    expect(statusLabel('interrupted')).toBe('已中断');
  });

  it('未知状态原样返回', () => {
    expect(statusLabel('weird')).toBe('weird');
  });
});

describe('stageLabel', () => {
  it('映射已知阶段', () => {
    expect(stageLabel('docking')).toBe('对接采样');
    expect(stageLabel('analyzing')).toBe('相互作用分析');
  });

  it('未知阶段原样返回', () => {
    expect(stageLabel('custom_stage')).toBe('custom_stage');
  });
});

describe('interactionTypeLabel', () => {
  it('映射已知相互作用类型', () => {
    expect(interactionTypeLabel('hbond')).toBe('氢键');
    expect(interactionTypeLabel('hydrophobic')).toBe('疏水作用');
    expect(interactionTypeLabel('pi_stacking')).toBe('π堆积');
    expect(interactionTypeLabel('salt_bridge')).toBe('盐桥');
  });

  it('未知类型原样返回', () => {
    expect(interactionTypeLabel('cation_pi')).toBe('cation_pi');
  });
});

describe('isTerminalStatus', () => {
  it('识别终态', () => {
    expect(isTerminalStatus('completed')).toBe(true);
    expect(isTerminalStatus('failed')).toBe(true);
    expect(isTerminalStatus('canceled')).toBe(true);
    expect(isTerminalStatus('interrupted')).toBe(true);
  });

  it('非终态', () => {
    expect(isTerminalStatus('queued')).toBe(false);
    expect(isTerminalStatus('running')).toBe(false);
  });
});

describe('parseResidue', () => {
  it('解析 残基名+序号:链 格式', () => {
    expect(parseResidue('ASP25:A')).toEqual({ resn: 'ASP', resi: 25, chain: 'A' });
  });

  it('无链时省略 chain', () => {
    expect(parseResidue('LEU102')).toEqual({ resn: 'LEU', resi: 102 });
  });

  it('负残基序号合法（PDB 允许）', () => {
    expect(parseResidue('ASP-25:B')).toEqual({ resn: 'ASP', resi: -25, chain: 'B' });
  });

  it('非法格式返回 null', () => {
    expect(parseResidue('')).toBeNull();
    expect(parseResidue('A25:')).toBeNull();
    expect(parseResidue('ASP25:')).toBeNull();
  });
});
