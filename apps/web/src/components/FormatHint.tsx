import { useState } from 'react';
import { useI18n } from '../i18n';

export type FormatHintKind = 'receptor' | 'ligand' | 'library';

const HINT_KEYS: Record<FormatHintKind, string[]> = {
  receptor: [
    'hint.receptor.format',
    'hint.receptor.cif',
    'hint.receptor.size',
    'hint.receptor.source',
  ],
  ligand: ['hint.ligand.format', 'hint.ligand.salt', 'hint.ligand.pubchem'],
  library: ['hint.library.smi', 'hint.library.csv', 'hint.library.source'],
};

const SMI_EXAMPLE = ['CC(=O)Oc1ccccc1C(=O)O aspirin', 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C caffeine', 'CCO ethanol'].join('\n');

interface FormatHintProps {
  kind: FormatHintKind;
}

export default function FormatHint({ kind }: FormatHintProps) {
  const { t } = useI18n();
  const [open, setOpen] = useState(false);

  return (
    <div className={open ? 'format-hint open' : 'format-hint'}>
      <button
        type="button"
        className="format-hint-toggle"
        aria-expanded={open}
        onClick={() => setOpen((v) => !v)}
      >
        <span className="format-hint-arrow">{open ? '▾' : '▸'}</span>
        {t('hint.toggle')}
      </button>
      {open && (
        <ul className="format-hint-body">
          {HINT_KEYS[kind].map((key) => (
            <li key={key}>{t(key)}</li>
          ))}
          {kind === 'library' && <pre className="format-hint-example">{SMI_EXAMPLE}</pre>}
        </ul>
      )}
    </div>
  );
}
