import * as d3ScaleChromatic from 'd3-scale-chromatic';
import * as d3Scale from 'd3-scale';
import { ColorSpec } from 'lib/data-map/view-layers';
import { valueType } from 'lib/helpers';

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

function invertColorScale<T>(colorScale: (t: number, n: number) => T) {
  return (i: number, n: number) => colorScale(1 - i, n);
}

export const VECTOR_COLOR_MAPS = valueType<ColorSpec>()({
  damages: {
    scale: d3Scale.scaleSequential,
    scheme: invertColorScale(d3ScaleChromatic.interpolateInferno),
    range: [0, 1000000],
    empty: '#ccc',
  },
  population: {
    scale: d3Scale.scaleSequential,
    scheme: d3ScaleChromatic.interpolateInferno,
    range: [0, 10000],
    empty: '#ccc',
  },
  adaptationAvoided: {
    scale: d3Scale.scaleSequential,
    scheme: d3ScaleChromatic.interpolateBlues,
    range: [0, 1000000],
    empty: '#ccc',
  },
  adaptationCost: {
    scale: d3Scale.scaleSequential,
    scheme: d3ScaleChromatic.interpolateGreens,
    range: [0, 1000000000],
    empty: '#ccc',
  },
});
