import type { DockRequest } from './types';

function fmt(n: number): string {
  return String(n);
}

export function renderVinaCli(
  req: DockRequest,
  options: { receptorFile?: string; ligandFile?: string; flex?: boolean } = {},
): string {
  const parts = [
    'vina',
    '--receptor',
    options.receptorFile ?? '<receptor.pdbqt>',
    '--ligand',
    options.ligandFile ?? '<ligand.pdbqt>',
    '--center_x',
    fmt(req.box.center[0]),
    '--center_y',
    fmt(req.box.center[1]),
    '--center_z',
    fmt(req.box.center[2]),
    '--size_x',
    fmt(req.box.size[0]),
    '--size_y',
    fmt(req.box.size[1]),
    '--size_z',
    fmt(req.box.size[2]),
    '--exhaustiveness',
    fmt(req.params.exhaustiveness),
    '--num_modes',
    fmt(req.params.num_modes),
    '--energy_range',
    fmt(req.params.energy_range),
    '--cpu',
    fmt(req.params.cpu),
    '--scoring',
    req.params.scoring,
  ];
  if (req.params.seed !== null) {
    parts.push('--seed', fmt(req.params.seed));
  }
  if (options.flex) {
    parts.push('--flex', '<flex.pdbqt>');
  }
  return parts.join(' ');
}
