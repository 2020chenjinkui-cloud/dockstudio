import { forwardRef, lazy, Suspense } from 'react';
import type { MolViewerHandle, MolViewerProps } from './MolViewer';

export type { MolViewerHandle, MolViewerProps, ResidueHighlight, ViewerPose } from './MolViewer';

const MolViewerImpl = lazy(() => import('./MolViewer'));

const MolViewerLazy = forwardRef<MolViewerHandle, MolViewerProps>(function MolViewerLazy(
  props,
  ref,
) {
  return (
    <Suspense fallback={<div className="mol-viewer mol-viewer-loading">3D…</div>}>
      <MolViewerImpl ref={ref} {...props} />
    </Suspense>
  );
});

export default MolViewerLazy;
