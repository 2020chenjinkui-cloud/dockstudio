import { useEffect, useRef } from 'react';
import * as $3Dmol from '3dmol';

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

interface MolViewerProps {
  receptorPdb: string | null;
  poses: ViewerPose[];
  selectedIndex: number | null;
  boxCenter: [number, number, number];
  boxSize: [number, number, number];
  highlightResidue?: ResidueHighlight | null;
}

export default function MolViewer({
  receptorPdb,
  poses,
  selectedIndex,
  boxCenter,
  boxSize,
  highlightResidue = null,
}: MolViewerProps) {
  const containerRef = useRef<HTMLDivElement>(null);
  const viewerRef = useRef<$3Dmol.GLViewer | null>(null);

  useEffect(() => {
    const container = containerRef.current;
    if (!container) return;
    if (!viewerRef.current) {
      viewerRef.current = $3Dmol.createViewer(container, {
        backgroundColor: '#151a23',
      });
    }
    const viewer = viewerRef.current;
    viewer.clear();

    const [cx, cy, cz] = boxCenter;
    const [sx, sy, sz] = boxSize;

    if (receptorPdb) {
      const receptor = viewer.addModel(receptorPdb, 'pdb');
      receptor.setStyle({}, { cartoon: { color: 'spectrum' } });
      receptor.setStyle(
        {
          hetflag: false,
          within: { distance: 5, sel: { x: cx, y: cy, z: cz } },
        },
        {
          cartoon: { color: 'spectrum' },
          stick: { colorscheme: 'whiteCarbon', radius: 0.15 },
        },
      );
      receptor.setStyle({ hetflag: true }, { stick: { colorscheme: 'greenCarbon', radius: 0.2 } });
      if (highlightResidue) {
        const sel: $3Dmol.AtomSelectionSpec = {
          resn: highlightResidue.resn,
          resi: highlightResidue.resi,
        };
        if (highlightResidue.chain) sel.chain = highlightResidue.chain;
        receptor.setStyle(sel, {
          cartoon: { color: 'spectrum' },
          stick: { colorscheme: 'yellowCarbon', radius: 0.3 },
          sphere: { colorscheme: 'yellowCarbon', scale: 0.35 },
        });
      }
    }

    if (receptorPdb || poses.some((p) => p.visible)) {
      viewer.addBox({
        center: { x: cx, y: cy, z: cz },
        dimensions: { w: sx, h: sy, d: sz },
        color: '#f2c14e',
        opacity: 0.3,
      });
    }

    for (const pose of poses) {
      if (!pose.visible) continue;
      const model = viewer.addModel(pose.sdf, 'sdf');
      if (pose.index === selectedIndex) {
        model.setStyle(
          {},
          {
            stick: { colorscheme: 'magentaCarbon', radius: 0.25 },
            sphere: { colorscheme: 'magentaCarbon', scale: 0.3 },
          },
        );
      } else {
        model.setStyle({}, { stick: { colorscheme: 'cyanCarbon', radius: 0.15 } });
      }
    }

    viewer.zoomTo();
    viewer.render();
  }, [receptorPdb, poses, selectedIndex, boxCenter, boxSize, highlightResidue]);

  useEffect(() => {
    const onResize = () => viewerRef.current?.resize();
    window.addEventListener('resize', onResize);
    return () => window.removeEventListener('resize', onResize);
  }, []);

  return <div ref={containerRef} className="mol-viewer" />;
}
