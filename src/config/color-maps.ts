import { scaleSequential } from 'd3-scale';
import { interpolateOrRd } from 'd3-scale-chromatic';

export const RASTER_COLOR_MAPS = {
  fluvial: {
    scheme: 'blues',
    range: [0, 10],
  },
  coastal: {
    scheme: 'greens',
    range: [0, 10],
  },
  surface: {
    scheme: 'purples',
    range: [0, 10],
  },
  cyclone: {
    scheme: 'reds',
    range: [0, 75],
  },
};

export const VECTOR_COLOR_MAPS = {
  damages: {
    scale: scaleSequential([0, 1000], interpolateOrRd),
    empty: '#ccc',
  },
};
