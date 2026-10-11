import type { DemoLibraryItem } from './types';

export const DEMO_PROJECT_NAME = '示例：1HSG';

export function demoLibraryToSmi(items: DemoLibraryItem[]): string {
  return items.map((item) => `${item.smiles} ${item.name}`.trimEnd()).join('\n');
}

export function demoLibraryToFile(name: string, items: DemoLibraryItem[]): File {
  return new File([demoLibraryToSmi(items)], `${name}.smi`, { type: 'text/plain' });
}
