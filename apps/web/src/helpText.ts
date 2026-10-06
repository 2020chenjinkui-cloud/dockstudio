export const PARAM_HELP = {
  exhaustiveness:
    '搜索彻底程度：值越大，构象搜索越充分、越耗时。默认 8 是精度与速度的平衡点；出最终结果前可提高到 16–32 复核。',
  num_modes:
    '最多输出多少个结合姿态。默认 9：足以覆盖主要结合模式，又不至于引入过多低分姿态。',
  cpu: '并行使用的 CPU 核数。默认 1；本机空闲时可调高以加速搜索。',
  scoring:
    '打分函数，用于评估姿态结合能力。默认 vina：AutoDock Vina 经典打分；引擎支持更多时可在下拉中切换。',
  boxCenter:
    '对接搜索空间的中心坐标（Å）。一般取共晶配体或已知口袋的质心；可用上一步共晶配体表的"定盒"自动填写。',
  boxSize:
    '对接搜索空间的边长（Å）。默认 22.5：覆盖常见口袋并留有余量；过大会稀释搜索、过小会漏掉结合姿态。',
} as const;

export const HINT_ADVANCED_MODE = 'M2 提供';
export const HINT_POCKET_SUGGESTION = 'M3 提供自动口袋预测';

export const DISCLAIMER =
  '计算结果仅供科研参考，对接打分不等于实测亲和力；任何结论需结合实验验证。';
