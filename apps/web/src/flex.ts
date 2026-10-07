export const FLEX_MAX = 10;

export function parseFlexInput(text: string): string[] {
  return text
    .split(/[\s,;，；\n]+/)
    .map((s) => s.trim())
    .filter(Boolean);
}

export interface FlexValidation {
  valid: string[];
  errors: string[];
}

export function validateFlexResidues(
  entries: string[],
  chains: string[],
  max: number = FLEX_MAX,
): FlexValidation {
  const valid: string[] = [];
  const errors: string[] = [];
  const seen = new Set<string>();
  for (const entry of entries) {
    const m = /^([A-Za-z0-9]+):(-?\d+)$/.exec(entry);
    if (!m) {
      errors.push(`${entry}：格式应为 链:残基号（如 A:42）`);
      continue;
    }
    const chain = m[1];
    if (chains.length > 0 && !chains.includes(chain)) {
      errors.push(`${entry}：受体中不存在链 ${chain}（可用：${chains.join(' / ')}）`);
      continue;
    }
    const normalized = `${chain}:${m[2]}`;
    if (seen.has(normalized)) {
      errors.push(`${entry}：重复`);
      continue;
    }
    seen.add(normalized);
    valid.push(normalized);
  }
  if (valid.length > max) {
    errors.push(`柔性残基最多 ${max} 个（当前 ${valid.length} 个）`);
  }
  return { valid, errors };
}
