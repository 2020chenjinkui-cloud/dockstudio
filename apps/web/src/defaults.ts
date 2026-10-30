export function defaultDockCpu(): number {
  const count =
    typeof navigator !== 'undefined' && navigator.hardwareConcurrency
      ? navigator.hardwareConcurrency
      : 4;
  return Math.max(2, Math.min(count, 8));
}
