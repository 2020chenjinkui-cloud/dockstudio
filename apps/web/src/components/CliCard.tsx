import { useMemo, useState } from 'react';
import { renderVinaCli } from '../cli';
import type { DockRequest } from '../types';

interface CliCardProps {
  request: DockRequest;
  flex?: boolean;
}

export default function CliCard({ request, flex = false }: CliCardProps) {
  const cmd = useMemo(() => renderVinaCli(request, { flex }), [request, flex]);
  const [copied, setCopied] = useState(false);

  const copy = async () => {
    try {
      await navigator.clipboard.writeText(cmd);
      setCopied(true);
      window.setTimeout(() => setCopied(false), 1500);
    } catch {
      // 剪贴板不可用时静默
    }
  };

  return (
    <div className="cli-card">
      <div className="cli-head">
        <h3 className="section-title">等效命令行</h3>
        <button onClick={copy}>{copied ? '已复制' : '复制'}</button>
      </div>
      <pre className="cli">
        <code>{cmd}</code>
      </pre>
    </div>
  );
}
