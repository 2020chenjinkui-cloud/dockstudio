import { describe, expect, it } from 'vitest';
import { renderVinaCli } from './cli';
import type { DockRequest } from './types';

const BASE: DockRequest = {
  engine: 'vina',
  receptor_id: 'r1',
  ligand_id: 'l1',
  box: { center: [1.5, -2, 30.25], size: [22.5, 22.5, 22.5] },
  params: {
    exhaustiveness: 8,
    num_modes: 9,
    energy_range: 3,
    cpu: 4,
    seed: null,
    scoring: 'vina',
  },
};

describe('renderVinaCli', () => {
  it('包含全部核心参数与占位文件名', () => {
    const cmd = renderVinaCli(BASE);
    expect(cmd).toContain('vina --receptor <receptor.pdbqt> --ligand <ligand.pdbqt>');
    expect(cmd).toContain('--center_x 1.5');
    expect(cmd).toContain('--center_y -2');
    expect(cmd).toContain('--center_z 30.25');
    expect(cmd).toContain('--size_x 22.5 --size_y 22.5 --size_z 22.5');
    expect(cmd).toContain('--exhaustiveness 8');
    expect(cmd).toContain('--num_modes 9');
    expect(cmd).toContain('--energy_range 3');
    expect(cmd).toContain('--cpu 4');
    expect(cmd).toContain('--scoring vina');
  });

  it('seed 为 null 时省略 --seed', () => {
    expect(renderVinaCli(BASE)).not.toContain('--seed');
  });

  it('seed 固定时带上 --seed', () => {
    const cmd = renderVinaCli({ ...BASE, params: { ...BASE.params, seed: 42 } });
    expect(cmd).toContain('--seed 42');
  });

  it('flex=true 时追加 --flex', () => {
    const cmd = renderVinaCli(BASE, { flex: true });
    expect(cmd).toContain('--flex <flex.pdbqt>');
    expect(renderVinaCli(BASE)).not.toContain('--flex');
  });
});
