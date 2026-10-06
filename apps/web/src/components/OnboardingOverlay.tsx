import { useState } from 'react';

const STEPS = [
  {
    title: '创建项目',
    body: '在顶部项目栏新建或切换项目。每个项目独立保存受体、配体与对接任务历史。',
  },
  {
    title: '按向导完成四步准备',
    body: '依次选择受体、配体、确定对接位点（可读取共晶配体自动定盒），然后提交运行。',
  },
  {
    title: '查看结果与报告',
    body: '结果页展示姿态打分、3D 相互作用与残基高亮，可一键下载自包含 HTML 报告。',
  },
];

export const ONBOARDING_KEY = 'dockstudio.onboarded.v1';

interface OnboardingOverlayProps {
  onDone: () => void;
}

export default function OnboardingOverlay({ onDone }: OnboardingOverlayProps) {
  const [index, setIndex] = useState(0);
  const step = STEPS[index];
  return (
    <div className="modal-overlay">
      <div className="modal onboarding">
        <div className="onboarding-dots">
          {STEPS.map((_, i) => (
            <span key={i} className={i === index ? 'dot-nav active' : 'dot-nav'} />
          ))}
        </div>
        <h3>{step.title}</h3>
        <p>{step.body}</p>
        <div className="modal-actions">
          <button onClick={onDone}>跳过</button>
          {index > 0 && <button onClick={() => setIndex(index - 1)}>上一步</button>}
          {index < STEPS.length - 1 ? (
            <button className="primary" onClick={() => setIndex(index + 1)}>
              下一步
            </button>
          ) : (
            <button className="primary" onClick={onDone}>
              开始使用
            </button>
          )}
        </div>
      </div>
    </div>
  );
}
