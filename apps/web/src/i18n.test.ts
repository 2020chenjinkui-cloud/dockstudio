import { describe, expect, it } from 'vitest';
import { DICTS, translate } from './i18n';

describe('i18n', () => {
  it('zh/en 两语言包键集合一致', () => {
    const zhKeys = Object.keys(DICTS.zh).sort();
    const enKeys = Object.keys(DICTS.en).sort();
    expect(enKeys).toEqual(zhKeys);
  });

  it('语言包没有空值', () => {
    for (const dict of Object.values(DICTS)) {
      for (const [key, value] of Object.entries(dict)) {
        expect(value, key).not.toBe('');
      }
    }
  });

  it('t() 正常取值', () => {
    expect(translate('zh', 'step.receptor')).toBe('受体');
    expect(translate('en', 'step.receptor')).toBe('Receptor');
  });

  it('t() 未知 key 优先 fallback，其次 key 本身', () => {
    expect(translate('en', 'status.mystery', 'mystery')).toBe('mystery');
    expect(translate('zh', 'status.mystery', 'raw_stage')).toBe('raw_stage');
    expect(translate('en', 'no.such.key')).toBe('no.such.key');
  });

  it('en 缺失某键时回退到 zh（运行时兜底）', () => {
    const dicts = DICTS as Record<string, Record<string, string>>;
    const saved = dicts.en['step.run'];
    delete dicts.en['step.run'];
    expect(translate('en', 'step.run')).toBe('运行');
    dicts.en['step.run'] = saved;
  });
});
