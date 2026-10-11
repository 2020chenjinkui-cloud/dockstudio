import { forwardRef, useEffect, useImperativeHandle, useMemo, useRef } from 'react';
import * as $3Dmol from '3dmol';
import { viewerContentKey } from './viewerLogic';

export interface ViewerPose {
  index: number;
  sdf: string;
  visible: boolean;
}

export interface ResidueHighlight {
  resn: string;
  resi: number;
  chain?: string;
}

export interface MolViewerHandle {
  capturePng: () => string | null;
}

export interface MolViewerProps {
  receptorPdb: string | null;
  poses: ViewerPose[];
  selectedIndex: number | null;
  boxCenter: [number, number, number];
  boxSize: [number, number, number];
  highlightResidue?: ResidueHighlight | null;
  crystalSdf?: string | null;
}

interface ModelMap {
  receptor: $3Dmol.GLModel | null;
  crystal: $3Dmol.GLModel | null;
  poses: Map<number, $3Dmol.GLModel>;
  poseSdfs: Map<number, string>;
}

function poseStyle(selected: boolean): $3Dmol.AtomStyleSpec {
  return selected
    ? {
        stick: { colorscheme: 'magentaCarbon', radius: 0.25 },
        sphere: { colorscheme: 'magentaCarbon', scale: 0.3 },
      }
    : { stick: { colorscheme: 'cyanCarbon', radius: 0.15 } };
}

function receptorStyles(
  receptor: $3Dmol.GLModel,
  boxCenter: [number, number, number],
  highlight: ResidueHighlight | null | undefined,
): void {
  const [cx, cy, cz] = boxCenter;
  receptor.setStyle({}, { cartoon: { color: 'spectrum' } });
  receptor.setStyle(
    { hetflag: false, within: { distance: 5, sel: { x: cx, y: cy, z: cz } } },
    { cartoon: { color: 'spectrum' }, stick: { colorscheme: 'whiteCarbon', radius: 0.15 } },
  );
  receptor.setStyle({ hetflag: true }, { stick: { colorscheme: 'greenCarbon', radius: 0.2 } });
  if (highlight) {
    const sel: $3Dmol.AtomSelectionSpec = { resn: highlight.resn, resi: highlight.resi };
    if (highlight.chain) sel.chain = highlight.chain;
    receptor.setStyle(sel, {
      cartoon: { color: 'spectrum' },
      stick: { colorscheme: 'yellowCarbon', radius: 0.3 },
      sphere: { colorscheme: 'yellowCarbon', scale: 0.35 },
    });
  }
}

const MolViewer = forwardRef<MolViewerHandle, MolViewerProps>(function MolViewer(props, ref) {
  const {
    receptorPdb,
    poses,
    selectedIndex,
    boxCenter,
    boxSize,
    highlightResidue = null,
    crystalSdf = null,
  } = props;

  const containerRef = useRef<HTMLDivElement>(null);
  const viewerRef = useRef<$3Dmol.GLViewer | null>(null);
  const modelsRef = useRef<ModelMap>({
    receptor: null,
    crystal: null,
    poses: new Map(),
    poseSdfs: new Map(),
  });
  const boxShapeRef = useRef<$3Dmol.GLShape | null>(null);
  const builtKeyRef = useRef<string>('');
  const hasZoomedRef = useRef(false);
  const propsRef = useRef(props);
  propsRef.current = props;

  const contentKey = useMemo(
    () => viewerContentKey(receptorPdb, poses, crystalSdf),
    [receptorPdb, poses, crystalSdf],
  );

  useImperativeHandle(
    ref,
    () => ({
      capturePng: () => viewerRef.current?.pngURI() ?? null,
    }),
    [],
  );

  useEffect(() => {
    const container = containerRef.current;
    if (!container) return;
    if (!viewerRef.current) {
      viewerRef.current = $3Dmol.createViewer(container, { backgroundColor: '#151a23' });
    }
    if (builtKeyRef.current === contentKey) return;
    builtKeyRef.current = contentKey;

    const viewer = viewerRef.current;
    const latest = propsRef.current;
    viewer.clear();
    const models = modelsRef.current;
    models.poses.clear();
    models.poseSdfs.clear();
    models.receptor = null;
    models.crystal = null;
    boxShapeRef.current = null;

    if (latest.receptorPdb) {
      models.receptor = viewer.addModel(latest.receptorPdb, 'pdb');
      receptorStyles(models.receptor, latest.boxCenter, latest.highlightResidue);
    }
    if (latest.crystalSdf) {
      models.crystal = viewer.addModel(latest.crystalSdf, 'sdf');
      models.crystal.setStyle(
        {},
        { stick: { color: '#d8d8d8', radius: 0.22, opacity: 0.55 } },
      );
    }
    for (const pose of latest.poses) {
      if (!pose.visible) continue;
      const model = viewer.addModel(pose.sdf, 'sdf');
      model.setStyle({}, poseStyle(pose.index === latest.selectedIndex));
      models.poses.set(pose.index, model);
      models.poseSdfs.set(pose.index, pose.sdf);
    }
    hasZoomedRef.current = false;
    viewer.zoomTo();
    hasZoomedRef.current = true;
    viewer.render();
  }, [contentKey]);

  useEffect(() => {
    const viewer = viewerRef.current;
    if (!viewer || builtKeyRef.current !== contentKey) return;
    const models = modelsRef.current;

    if (models.receptor) {
      receptorStyles(models.receptor, boxCenter, highlightResidue);
    }
    for (const pose of poses) {
      const existing = models.poses.get(pose.index);
      if (pose.visible && !existing) {
        const sdf = models.poseSdfs.get(pose.index) ?? pose.sdf;
        const model = viewer.addModel(sdf, 'sdf');
        models.poses.set(pose.index, model);
        models.poseSdfs.set(pose.index, sdf);
        model.setStyle({}, poseStyle(pose.index === selectedIndex));
      } else if (!pose.visible && existing) {
        viewer.removeModel(existing);
        models.poses.delete(pose.index);
      } else if (existing) {
        existing.setStyle({}, poseStyle(pose.index === selectedIndex));
      }
    }

    if (boxShapeRef.current) {
      viewer.removeShape(boxShapeRef.current);
      boxShapeRef.current = null;
    }
    const hasAnything =
      models.receptor !== null || models.poses.size > 0 || models.crystal !== null;
    if (hasAnything) {
      const [cx, cy, cz] = boxCenter;
      const [sx, sy, sz] = boxSize;
      boxShapeRef.current = viewer.addBox({
        center: { x: cx, y: cy, z: cz },
        dimensions: { w: sx, h: sy, d: sz },
        color: '#f2c14e',
        opacity: 0.3,
      });
    }
    viewer.render();
  }, [contentKey, poses, selectedIndex, boxCenter, boxSize, highlightResidue]);

  useEffect(() => {
    const onResize = () => viewerRef.current?.resize();
    window.addEventListener('resize', onResize);
    return () => window.removeEventListener('resize', onResize);
  }, []);

  return <div ref={containerRef} className="mol-viewer" />;
});

export default MolViewer;
