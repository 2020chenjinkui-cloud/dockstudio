import type { ComparisonRow } from './types';

export function toggleCompareId(ids: string[], id: string): string[] {
  return ids.includes(id) ? ids.filter((x) => x !== id) : [...ids, id];
}

export function bestAffinity(rows: ComparisonRow[]): number | null {
  const values = rows.map((r) => r.top_affinity).filter((v): v is number => v !== null);
  if (values.length === 0) return null;
  return Math.min(...values);
}

function sanitizeCsvCell(value: string): string {
  // Prevent spreadsheet formula/DDE injection when the export is opened in Excel.
  return /^[=+\-@\t\r]/.test(value) ? `'${value}` : value;
}

function escCsv(value: string): string {
  const safe = sanitizeCsvCell(value);
  return /[",\n]/.test(safe) ? `"${safe.replace(/"/g, '""')}"` : safe;
}

export function toCompareCsv(rows: ComparisonRow[]): string {
  const header =
    'job_id,ligand_name,smiles,engine,scoring,exhaustiveness,top_affinity,n_poses,status';
  const lines = rows.map((r) =>
    [
      r.job_id,
      escCsv(r.ligand_name ?? ''),
      escCsv(r.smiles ?? ''),
      r.engine,
      r.scoring,
      String(r.exhaustiveness),
      r.top_affinity !== null ? r.top_affinity.toFixed(2) : '',
      String(r.n_poses),
      r.status,
    ].join(','),
  );
  return [header, ...lines].join('\n');
}
