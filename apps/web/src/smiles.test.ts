import { describe, expect, it } from 'vitest';
import { checkSmilesQuick } from './smiles';

describe('checkSmilesQuick', () => {
  it('空输入不出卡', () => {
    expect(checkSmilesQuick('')).toBeNull();
    expect(checkSmilesQuick('   ')).toBeNull();
  });

  it('常规 SMILES 判定为合理', () => {
    const r = checkSmilesQuick('CC(=O)Oc1ccccc1C(=O)O');
    expect(r).toEqual({ verdict: 'ok', messageKey: 'smiles.ok' });
  });

  it('纯数字拒绝', () => {
    const r = checkSmilesQuick('12345');
    expect(r?.verdict).toBe('error');
    expect(r?.messageKey).toBe('smiles.numeric');
  });

  it('含空格提示', () => {
    const r = checkSmilesQuick('CCO ethanol');
    expect(r?.verdict).toBe('warn');
    expect(r?.messageKey).toBe('smiles.space');
  });

  it('超长提示', () => {
    const r = checkSmilesQuick('C'.repeat(1001));
    expect(r?.verdict).toBe('warn');
    expect(r?.messageKey).toBe('smiles.tooLong');
  });

  it('非法字符拒绝', () => {
    const r = checkSmilesQuick('CCO!bad');
    expect(r?.verdict).toBe('error');
    expect(r?.messageKey).toBe('smiles.invalid');
  });

  it('盐分隔符提醒拆分', () => {
    const r = checkSmilesQuick('c1ccccc1.[Na+]');
    expect(r?.verdict).toBe('warn');
    expect(r?.messageKey).toBe('smiles.salt');
  });
});
