export const PARAM_HELP = {
  exhaustiveness:
    '搜索彻底程度：值越大，构象搜索越充分、越耗时。默认 8 是精度与速度的平衡点；出最终结果前可提高到 16–32 复核。',
  num_modes:
    '最多输出多少个结合姿态。默认 9：足以覆盖主要结合模式，又不至于引入过多低分姿态。',
  cpu: '并行使用的 CPU 核数。默认 1；本机空闲时可调高以加速搜索。',
  scoring:
    '打分函数，用于评估姿态结合能力。默认 vina：AutoDock Vina 经典打分；引擎支持更多时可在下拉中切换。',
  energyRange:
    '输出姿态的能量范围（kcal/mol）：只保留与最佳姿态能量差在此范围内的姿态。默认 3.0。',
  seed: '随机种子：留空表示随机；填入固定整数可复现实验结果。',
  flexResidues:
    '柔性残基：对接中允许侧链自由运动的残基。格式 链:残基号（如 A:42），多个用空格或逗号分隔，最多 10 个；留空则刚性对接。',
  boxCenter:
    '对接搜索空间的中心坐标（Å）。一般取共晶配体或已知口袋的质心；可用上一步共晶配体表的"定盒"自动填写。',
  boxSize:
    '对接搜索空间的边长（Å）。默认 22.5：覆盖常见口袋并留有余量；过大会稀释搜索、过小会漏掉结合姿态。',
  batchLibrary:
    '配体库文件：.sdf（多分子）/ .smi（每行 SMILES [名称]）/ .csv（含 smiles 列）。单条解析失败不影响其余条目，上限 50000 条。',
} as const;

export const DISCLAIMER =
  '计算结果仅供科研参考，对接打分不等于实测亲和力；任何结论需结合实验验证。';

export const AI_HELP = {
  what:
    'AI 复核用 Boltz-2 对选定姿态独立预测 log(IC50) 与结合概率，与 vina 打分互不依赖；数值不一致属正常，用于交叉验证打分可靠性。',
  logIC50:
    'log10(IC50, μM)：AI 预测的半数抑制浓度对数值。越低预测亲和力越强；与 vina affinity 是两套独立模型，不保证一致。',
  binderProb:
    'binder_probability：AI 判断该小分子为真实结合者的概率（0–100%）。越高越可信，但仍是计算预测，需实验验证。',
} as const;

export const AI_DOC_HINT = '安装与配置指引见仓库 docs/ai.md。';
