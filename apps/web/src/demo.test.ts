import { describe, expect, it } from 'vitest';
import { demoLibraryToSmi } from './demo';

describe('demoLibraryToSmi', () => {
  it('每行 SMILES 空格 名称', () => {
    const smi = demoLibraryToSmi([
      { name: 'aspirin', smiles: 'CC(=O)Oc1ccccc1C(=O)O' },
      { name: 'ethanol', smiles: 'CCO' },
    ]);
    expect(smi.split('\n')).toEqual([
      'CC(=O)Oc1ccccc1C(=O)O aspirin',
      'CCO ethanol',
    ]);
  });

  it('空库输出空字符串', () => {
    expect(demoLibraryToSmi([])).toBe('');
  });

  it('名称为空时行尾不带空格', () => {
    expect(demoLibraryToSmi([{ name: '', smiles: 'CCO' }])).toBe('CCO');
  });
});
