import { describe, expect, it } from 'vitest';
import { crystalBadge, formatCrystalRmsd, viewerContentKey } from './viewerLogic';

const poseA = { index: 1, sdf: 'a'.repeat(100) };
const poseB = { index: 2, sdf: 'b'.repeat(200) };

describe('viewerContentKey', () => {
  it('同内容换数组实例 key 不变（增量路径保留相机）', () => {
    const k1 = viewerContentKey('PDBTEXT', [poseA, poseB]);
    const k2 = viewerContentKey('PDBTEXT', [
      { index: 2, sdf: 'b'.repeat(200) },
      { index: 1, sdf: 'a'.repeat(100) },
    ]);
    expect(k1).toBe(k2);
  });

  it('受体或姿态内容变化时 key 变化（触发全量重建）', () => {
    const base = viewerContentKey('PDBTEXT', [poseA]);
    expect(viewerContentKey('PDBTEXT2', [poseA])).not.toBe(base);
    expect(viewerContentKey('PDBTEXT', [poseA, poseB])).not.toBe(base);
    expect(viewerContentKey('PDBTEXT', [{ index: 1, sdf: 'x'.repeat(50) }])).not.toBe(base);
  });

  it('可见性/选中不参与 key（只走增量样式）', () => {
    const withVisible = [{ index: 1, sdf: 'a'.repeat(100), visible: false }];
    const plain = [{ index: 1, sdf: 'a'.repeat(100), visible: true }];
    expect(viewerContentKey(null, withVisible)).toBe(viewerContentKey(null, plain));
  });

  it('共晶配体进入 key', () => {
    const without = viewerContentKey('P', [poseA], null);
    const withCrystal = viewerContentKey('P', [poseA], 'SDFCRYSTAL');
    expect(withCrystal).not.toBe(without);
  });
});

describe('crystalBadge / formatCrystalRmsd', () => {
  it('null → none（列显示 —）', () => {
    expect(crystalBadge(null)).toBe('none');
    expect(formatCrystalRmsd(undefined)).toBeNull();
    expect(formatCrystalRmsd({ pose_index: 1, rmsd: null, reference: null })).toBeNull();
  });

  it('< 2 Å 判定通过', () => {
    expect(crystalBadge(1.85)).toBe('pass');
    expect(crystalBadge(2)).toBe('plain');
    expect(crystalBadge(3.4)).toBe('plain');
  });
});
