import { useI18n } from '../i18n';
import { crystalBadge } from '../viewerLogic';

interface CrystalRmsdCellProps {
  rmsd: number | null;
}

export default function CrystalRmsdCell({ rmsd }: CrystalRmsdCellProps) {
  const { t } = useI18n();
  const badge = crystalBadge(rmsd);
  if (badge === 'none') {
    return <span className="muted">—</span>;
  }
  return (
    <span className={badge === 'pass' ? 'crystal-pass' : undefined}>
      {rmsd!.toFixed(2)}
      {badge === 'pass' && <span className="badge qc-ok crystal-badge">{t('crystal.pass')}</span>}
    </span>
  );
}
