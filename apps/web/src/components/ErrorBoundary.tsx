import { Component, type ReactNode } from 'react';
import { translate, LANG_KEY, type Lang } from '../i18n';

interface ErrorBoundaryState {
  error: Error | null;
}

export default class ErrorBoundary extends Component<{ children: ReactNode }, ErrorBoundaryState> {
  constructor(props: { children: ReactNode }) {
    super(props);
    this.state = { error: null };
  }

  static getDerivedStateFromError(error: Error): ErrorBoundaryState {
    return { error };
  }

  componentDidCatch(error: Error, info: unknown): void {
    console.error('[ErrorBoundary]', error, info);
  }

  render() {
    if (this.state.error) {
      const lang: Lang = localStorage.getItem(LANG_KEY) === 'en' ? 'en' : 'zh';
      return (
        <div className="gate">
          <div className="card gate-card">
            <h2>{translate(lang, 'err.title')}</h2>
            <p className="muted">{translate(lang, 'err.desc')}</p>
            <pre className="log err-detail">
              {this.state.error.message}
            </pre>
            <div className="modal-actions">
              <button className="primary" onClick={() => window.location.reload()}>
                {translate(lang, 'err.reload')}
              </button>
            </div>
          </div>
        </div>
      );
    }
    return this.props.children;
  }
}
