import * as d3 from 'd3-scale-chromatic';

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

function invertColorScale(colorScale) {
  return (i, n) => colorScale(1 - i, n);
}

export const VECTOR_COLOR_MAPS = {
  damages: {
    scale: invertColorScale(d3.interpolateInferno),
    range: [0, 1000000],
    empty: '#ccc',
  },
};
