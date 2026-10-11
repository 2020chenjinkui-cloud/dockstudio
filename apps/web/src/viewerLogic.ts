import type { CrystalRmsd } from './types';

export interface ViewerPoseLike {
  index: number;
  sdf: string;
}

export function viewerContentKey(
  receptorPdb: string | null,
  poses: ViewerPoseLike[],
  crystalSdf: string | null = null,
): string {
  const poseSig = poses
    .map((p) => `${p.index}:${p.sdf.length}`)
    .sort()
    .join(',');
  return [
    receptorPdb === null ? 'none' : `pdb:${receptorPdb.length}`,
    poseSig === '' ? 'nopose' : `poses:${poseSig}`,
    crystalSdf === null ? 'nocrystal' : `crystal:${crystalSdf.length}`,
  ].join('|');
}

export type CrystalBadge = 'pass' | 'plain' | 'none';

export function crystalBadge(rmsd: number | null): CrystalBadge {
  if (rmsd === null) return 'none';
  return rmsd < 2 ? 'pass' : 'plain';
}

export function formatCrystalRmsd(entry: CrystalRmsd | undefined): number | null {
  return entry?.rmsd ?? null;
}
