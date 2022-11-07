import { d3Scale, d3ScaleChromatic, invertColorScale } from '@/lib/data-map/color-maps';
import { ColorSpec } from '@/lib/data-map/view-layers';

export const DAMAGE_COLORMAP: ColorSpec = {
  scale: d3Scale.scaleSequential,
  scheme: invertColorScale(d3ScaleChromatic.interpolateInferno),
  range: [0, 1e6],
  empty: '#ccc',
};
