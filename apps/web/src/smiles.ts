export type SmilesQuickVerdict = 'ok' | 'warn' | 'error';

export interface SmilesQuickCheck {
  verdict: SmilesQuickVerdict;
  messageKey: string;
}

const SMILES_CHARSET = /^[A-Za-z0-9()\[\]@+\-=#$:/\.%*]+$/;

export function checkSmilesQuick(input: string): SmilesQuickCheck | null {
  const s = input.trim();
  if (s === '') return null;
  if (/^\d+$/.test(s)) {
    return { verdict: 'error', messageKey: 'smiles.numeric' };
  }
  if (/\s/.test(s)) {
    return { verdict: 'warn', messageKey: 'smiles.space' };
  }
  if (s.length > 1000) {
    return { verdict: 'warn', messageKey: 'smiles.tooLong' };
  }
  if (!SMILES_CHARSET.test(s)) {
    return { verdict: 'error', messageKey: 'smiles.invalid' };
  }
  if (s.includes('.')) {
    return { verdict: 'warn', messageKey: 'smiles.salt' };
  }
  return { verdict: 'ok', messageKey: 'smiles.ok' };
}
