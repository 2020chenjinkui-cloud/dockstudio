import { useI18n } from '../../i18n';
import type { DockParams, EngineManifest } from '../../types';
import Tooltip from '../Tooltip';

interface ParamsPanelProps {
  params: DockParams;
  engines: EngineManifest[];
  engine: string;
  advanced: boolean;
  onEngineChange: (engine: string) => void;
  onParamsChange: (params: DockParams) => void;
}

function toNumber(value: string): number {
  const parsed = Number.parseFloat(value);
  return Number.isNaN(parsed) ? 0 : parsed;
}

export default function ParamsPanel({
  params,
  engines,
  engine,
  advanced,
  onEngineChange,
  onParamsChange,
}: ParamsPanelProps) {
  const { t } = useI18n();
  const selected = engines.find((item) => item.name === engine);
  const scoringOptions =
    selected && selected.scoring_functions.length > 0
      ? selected.scoring_functions
      : [params.scoring];
  const engineOptions = engines.length > 0 ? engines : [];

  const chooseEngine = (name: string) => {
    const next = engines.find((item) => item.name === name);
    const scoring =
      next && next.scoring_functions.includes(params.scoring)
        ? params.scoring
        : (next?.scoring_functions[0] ?? params.scoring);
    onParamsChange({ ...params, scoring });
    onEngineChange(name);
  };

  return (
    <div className="grid-form">
      <label>
        <span className="control-row">
          引擎
          <Tooltip text={t('help.scoring')} />
        </span>
        <select value={engine} onChange={(event) => chooseEngine(event.target.value)}>
          {engineOptions.length === 0 ? (
            <option value={engine}>{engine}</option>
          ) : (
            engineOptions.map((item) => (
              <option key={item.name} value={item.name} disabled={!item.available}>
                {item.name}
                {item.available
                  ? item.version && item.version !== 'unknown'
                    ? ` · ${item.version}`
                    : ''
                  : '（不可用）'}
              </option>
            ))
          )}
        </select>
      </label>
      <label>
        <span className="control-row">
          scoring
          <Tooltip text={t('help.scoring')} />
        </span>
        <select
          value={params.scoring}
          onChange={(event) => onParamsChange({ ...params, scoring: event.target.value })}
        >
          {scoringOptions.map((item) => (
            <option key={item} value={item}>
              {item}
            </option>
          ))}
        </select>
      </label>
      <label>
        <span className="control-row">
          exhaustiveness (1-64)
          <Tooltip text={t('help.exhaustiveness')} />
        </span>
        <input
          type="number"
          min={1}
          max={64}
          value={params.exhaustiveness}
          onChange={(event) =>
            onParamsChange({
              ...params,
              exhaustiveness: Math.min(64, Math.max(1, toNumber(event.target.value))),
            })
          }
        />
      </label>
      <label>
        <span className="control-row">
          num_modes (1-100)
          <Tooltip text={t('help.num_modes')} />
        </span>
        <input
          type="number"
          min={1}
          max={100}
          value={params.num_modes}
          onChange={(event) =>
            onParamsChange({
              ...params,
              num_modes: Math.min(100, Math.max(1, toNumber(event.target.value))),
            })
          }
        />
      </label>
      <label>
        <span className="control-row">
          cpu
          <Tooltip text={t('help.cpu')} />
        </span>
        <input
          type="number"
          min={1}
          value={params.cpu}
          onChange={(event) =>
            onParamsChange({ ...params, cpu: Math.max(1, toNumber(event.target.value)) })
          }
        />
      </label>
      {advanced && (
        <>
          <label>
            <span className="control-row">
              energy_range
              <Tooltip text={t('help.energyRange')} />
            </span>
            <input
              type="number"
              step={0.5}
              min={0.1}
              value={params.energy_range}
              onChange={(event) =>
                onParamsChange({ ...params, energy_range: toNumber(event.target.value) })
              }
            />
          </label>
          <label>
            <span className="control-row">
              seed（留空随机）
              <Tooltip text={t('help.seed')} />
            </span>
            <input
              type="text"
              placeholder="随机"
              value={params.seed === null ? '' : String(params.seed)}
              onChange={(event) => {
                const value = event.target.value.trim();
                onParamsChange({
                  ...params,
                  seed: value === '' ? null : Math.trunc(toNumber(value)),
                });
              }}
            />
          </label>
        </>
      )}
    </div>
  );
}
