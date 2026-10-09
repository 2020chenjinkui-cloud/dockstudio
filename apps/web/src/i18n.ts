import { createContext, useContext } from 'react';
import { PARAM_HELP } from './helpText';

export type Lang = 'zh' | 'en';

export const LANG_KEY = 'dockstudio.lang';

const zh = {
  // 向导
  'step.receptor': '受体',
  'step.ligand': '配体',
  'step.ligandLibrary': '配体库',
  'step.site': '对接位点',
  'step.run': '运行',
  'step.results': '结果',
  'nav.prev': '上一步',
  'nav.next': '下一步',
  'mode.single': '单一对接',
  'mode.batch': '批量筛选',
  'adv.toggle': '高级模式',
  'adv.on': '高级模式：开',
  // 状态徽章
  'status.queued': '排队中',
  'status.running': '运行中',
  'status.completed': '已完成',
  'status.failed': '失败',
  'status.canceled': '已取消',
  'status.interrupted': '已中断',
  // 阶段映射
  'stage.queued': '排队等待',
  'stage.preparing': '结构准备',
  'stage.docking': '对接采样',
  'stage.scoring': '打分排序',
  'stage.exporting': '导出姿态',
  'stage.analyzing': '相互作用分析',
  'stage.finished': '完成',
  'stage.done': '完成',
  // 项目栏
  'bar.new': '新建项目',
  'bar.history': '任务历史',
  'bar.delete': '删除项目',
  'bar.create': '创建',
  'bar.cancel': '取消',
  'bar.close': '关闭',
  'bar.backendOk': '后端已连接',
  'bar.backendDown': '后端未启动',
  'bar.namePlaceholder': '新项目名称',
  'bar.none': '未选择项目',
  'bar.jobsEmpty': '本项目还没有任务。',
  // 统计
  'stats.projects': '项目',
  'stats.completed': '已完成任务',
  'stats.avg': '平均对接',
  // 通用
  'common.loading': '加载中…',
  'common.refresh': '刷新',
  // 姿态质检
  'qc.ok': '合格',
  'qc.issue': '有问题',
  'qc.onlyOk': '只看合格',
  'qc.popoverTitle': '姿态质检',
  'qc.severity.info': '提示',
  'qc.severity.warning': '警告',
  'qc.severity.error': '严重',
  // 相互作用 / 2D 图
  'tab.interactions': '相互作用',
  'tab.diagram2d': '2D 示意图',
  'diagram.error': '2D 示意图加载失败（该姿态可能暂无示意图）',
  // 报告
  'report.download': '下载报告（HTML）',
  'report.dialogTitle': '下载报告',
  'report.sections': '包含章节',
  'report.customTitle': '自定义标题（可选，不超过 80 字）',
  'report.confirm': '下载',
  'report.section.snapshot': '参数快照',
  'report.section.scores': '打分表',
  'report.section.interactions': '相互作用',
  'report.section.diagram': '2D 示意图',
  'report.section.viewer3d': '3D 视图',
  'report.section.disclaimer': '免责声明',
  // 对比
  'cmp.button': '对比',
  'cmp.title': '结果对比',
  'cmp.hint': '在历史抽屉中勾选 2 个及以上已完成任务后可对比',
  'cmp.col.ligand': '配体',
  'cmp.col.smiles': 'SMILES',
  'cmp.col.engine': '引擎',
  'cmp.col.scoring': '打分函数',
  'cmp.col.exhaustiveness': 'exhaustiveness',
  'cmp.col.affinity': 'top affinity',
  'cmp.col.poses': '姿态数',
  'cmp.col.status': '状态',
  'cmp.export': '导出 CSV',
  'cmp.best': '最优',
  // AI 复核
  'ai.button': 'AI 复核',
  'ai.title': 'AI 复核（Boltz-2）',
  'ai.what': 'Boltz-2 独立预测 log(IC50) 与结合概率，用于交叉验证 vina 打分。',
  'ai.poses': '复核姿态（默认 top 3）',
  'ai.msa': '使用公共 MSA 服务器（ColabFold，更快）',
  'ai.msaWarn': '蛋白序列将发送至公共 ColabFold MSA 服务器',
  'ai.msaWarnHistory': '本次复核的蛋白序列已发送至公共 MSA 服务器',
  'ai.unavailableDocs': '安装与配置指引见仓库 docs/ai.md。',
  'ai.submit': '提交复核',
  'ai.submitting': '提交中…',
  'ai.running': '复核中',
  'ai.cancel': '取消复核',
  'ai.history': '复核记录',
  'ai.noError': '至少需要勾选一个姿态',
  'ai.posesCount': '姿态数',
  'ai.th.logIC50': 'logIC50',
  'ai.th.binder': 'binder%',
  'ai.note.logIC50': 'log10(IC50 μM)：AI 独立预测，与 vina 分数不保证一致',
  'ai.note.binder': 'binder 概率：AI 判断为真实结合者的置信度，AI 独立预测',
  'ai.mixNotice': 'AI 复核为独立预测，与 vina affinity 并排仅供参考',
  // 帮助 tooltip
  'help.exhaustiveness': PARAM_HELP.exhaustiveness,
  'help.num_modes': PARAM_HELP.num_modes,
  'help.cpu': PARAM_HELP.cpu,
  'help.scoring': PARAM_HELP.scoring,
  'help.energyRange': PARAM_HELP.energyRange,
  'help.seed': PARAM_HELP.seed,
  'help.flexResidues': PARAM_HELP.flexResidues,
  'help.boxCenter': PARAM_HELP.boxCenter,
  'help.boxSize': PARAM_HELP.boxSize,
  'help.batchLibrary': PARAM_HELP.batchLibrary,
} as const;

export type I18nKey = keyof typeof zh;

const en: Record<I18nKey, string> = {
  'step.receptor': 'Receptor',
  'step.ligand': 'Ligand',
  'step.ligandLibrary': 'Ligand library',
  'step.site': 'Binding site',
  'step.run': 'Run',
  'step.results': 'Results',
  'nav.prev': 'Back',
  'nav.next': 'Next',
  'mode.single': 'Single dock',
  'mode.batch': 'Batch screen',
  'adv.toggle': 'Advanced',
  'adv.on': 'Advanced: ON',
  'status.queued': 'Queued',
  'status.running': 'Running',
  'status.completed': 'Completed',
  'status.failed': 'Failed',
  'status.canceled': 'Canceled',
  'status.interrupted': 'Interrupted',
  'stage.queued': 'Queued',
  'stage.preparing': 'Preparing',
  'stage.docking': 'Docking',
  'stage.scoring': 'Scoring',
  'stage.exporting': 'Exporting poses',
  'stage.analyzing': 'Analyzing interactions',
  'stage.finished': 'Done',
  'stage.done': 'Done',
  'bar.new': 'New project',
  'bar.history': 'History',
  'bar.delete': 'Delete project',
  'bar.create': 'Create',
  'bar.cancel': 'Cancel',
  'bar.close': 'Close',
  'bar.backendOk': 'Backend connected',
  'bar.backendDown': 'Backend offline',
  'bar.namePlaceholder': 'New project name',
  'bar.none': 'No project selected',
  'bar.jobsEmpty': 'No jobs in this project yet.',
  'stats.projects': 'projects',
  'stats.completed': 'completed jobs',
  'stats.avg': 'avg dock',
  'common.loading': 'Loading…',
  'common.refresh': 'Refresh',
  'qc.ok': 'Pass',
  'qc.issue': 'Issues',
  'qc.onlyOk': 'Pass only',
  'qc.popoverTitle': 'Pose check',
  'qc.severity.info': 'Info',
  'qc.severity.warning': 'Warning',
  'qc.severity.error': 'Error',
  'tab.interactions': 'Interactions',
  'tab.diagram2d': '2D diagram',
  'diagram.error': 'Failed to load the 2D diagram (it may not exist for this pose)',
  'report.download': 'Download report (HTML)',
  'report.dialogTitle': 'Download report',
  'report.sections': 'Sections',
  'report.customTitle': 'Custom title (optional, up to 80 chars)',
  'report.confirm': 'Download',
  'report.section.snapshot': 'Snapshot',
  'report.section.scores': 'Scores',
  'report.section.interactions': 'Interactions',
  'report.section.diagram': '2D diagram',
  'report.section.viewer3d': '3D viewer',
  'report.section.disclaimer': 'Disclaimer',
  'cmp.button': 'Compare',
  'cmp.title': 'Compare results',
  'cmp.hint': 'Tick 2+ completed jobs in the history drawer to compare',
  'cmp.col.ligand': 'Ligand',
  'cmp.col.smiles': 'SMILES',
  'cmp.col.engine': 'Engine',
  'cmp.col.scoring': 'Scoring',
  'cmp.col.exhaustiveness': 'exhaustiveness',
  'cmp.col.affinity': 'top affinity',
  'cmp.col.poses': 'Poses',
  'cmp.col.status': 'Status',
  'cmp.export': 'Export CSV',
  'cmp.best': 'best',
  'ai.button': 'AI rescore',
  'ai.title': 'AI rescore (Boltz-2)',
  'ai.what': 'Boltz-2 independently predicts log(IC50) and binder probability to cross-validate vina scores.',
  'ai.poses': 'Poses to rescore (top 3 by default)',
  'ai.msa': 'Use the public MSA server (ColabFold, faster)',
  'ai.msaWarn': 'The protein sequence will be sent to the public ColabFold MSA server',
  'ai.msaWarnHistory': 'The protein sequence was sent to the public MSA server for this rescore',
  'ai.unavailableDocs': 'See docs/ai.md in the repository for setup instructions.',
  'ai.submit': 'Start rescore',
  'ai.submitting': 'Submitting…',
  'ai.running': 'Rescoring',
  'ai.cancel': 'Cancel rescore',
  'ai.history': 'Rescore history',
  'ai.noError': 'Select at least one pose',
  'ai.posesCount': 'poses',
  'ai.th.logIC50': 'logIC50',
  'ai.th.binder': 'binder%',
  'ai.note.logIC50': 'log10(IC50 μM): independent AI prediction, may disagree with vina scores',
  'ai.note.binder': 'Binder probability: confidence that the ligand truly binds; independent AI prediction',
  'ai.mixNotice': 'AI rescore is an independent prediction shown next to vina affinity for reference only',
  'help.exhaustiveness':
    'Thoroughness of conformational search: higher values are more exhaustive but slower. Default 8 balances accuracy and speed; raise to 16–32 for final confirmations.',
  'help.num_modes':
    'Maximum number of binding poses to output. Default 9 covers major binding modes without too many low-score poses.',
  'help.cpu': 'CPU cores used in parallel. Default 1; raise it when the machine is idle to speed up the search.',
  'help.scoring':
    'Scoring function used to evaluate poses. Default vina, the classic AutoDock Vina score; switch in the dropdown when the engine provides more.',
  'help.energyRange':
    'Energy range of reported poses (kcal/mol): only poses within this distance from the best pose are kept. Default 3.0.',
  'help.seed': 'Random seed: leave empty for random; set a fixed integer for reproducible runs.',
  'help.flexResidues':
    'Flexible residues whose side chains may move during docking. Format chain:resno (e.g. A:42), separated by spaces or commas, at most 10; leave empty for rigid docking.',
  'help.boxCenter':
    'Center of the docking search box (Å). Usually the centroid of a cocrystal ligand or a known pocket; use the cocrystal table in the previous step to auto-fill.',
  'help.boxSize':
    'Edge length of the docking search box (Å). Default 22.5 covers typical pockets with margin; too large dilutes sampling, too small misses poses.',
  'help.batchLibrary':
    'Ligand library file: .sdf (multi-molecule) / .smi (SMILES [name] per line) / .csv (with a smiles column). A failed entry does not stop the rest; up to 50000 entries.',
};

export const DICTS: Record<Lang, Record<I18nKey, string>> = { zh, en };

export function translate(lang: Lang, key: string, fallback?: string): string {
  const dict: Record<string, string> = DICTS[lang];
  return dict[key] ?? DICTS.zh[key as I18nKey] ?? fallback ?? key;
}

export interface I18nValue {
  lang: Lang;
  setLang: (lang: Lang) => void;
  t: (key: string, fallback?: string) => string;
}

export const I18nContext = createContext<I18nValue>({
  lang: 'zh',
  setLang: () => undefined,
  t: (key, fallback) => translate('zh', key, fallback),
});

export function useI18n(): I18nValue {
  return useContext(I18nContext);
}
