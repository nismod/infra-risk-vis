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

function invertColorScale<T>(colorScale: (t: number) => T) {
  return (i: number) => colorScale(1 - i);
}

function discardSides<T>(interpolator: (t: number) => T, cutStart: number, cutEnd: number = 0) {
  return (i: number) => {
    const t = i * (1 - cutStart - cutEnd) + cutStart;
    return interpolator(t);
  };
}

export const VECTOR_COLOR_MAPS = valueType<ColorSpec>()({
  damages: {
    scale: d3Scale.scaleSequential,
    scheme: invertColorScale(d3ScaleChromatic.interpolateInferno),
    range: [0, 1e6],
    empty: '#ccc',
  },
  population: {
    scale: d3Scale.scaleSequential,
    scheme: invertColorScale(d3ScaleChromatic.interpolateInferno),
    range: [0, 1e4],
    empty: '#ccc',
  },
  adaptationAvoided: {
    scale: d3Scale.scaleSequential,
    scheme: discardSides(d3ScaleChromatic.interpolateBlues, 0.2, 0.2),
    range: [0, 1e7],
    empty: '#ccc',
  },
  adaptationCost: {
    scale: d3Scale.scaleSequential,
    scheme: discardSides(d3ScaleChromatic.interpolateGreens, 0.2, 0.2),
    range: [0, 1e9],
    empty: '#ccc',
  },
  costBenefitRatio: {
    scale: d3Scale.scaleSequential,
    scheme: invertColorScale(d3ScaleChromatic.interpolateViridis),
    range: [1, 10],
    empty: '#ccc',
  },
  terrestrialElevation: {
    scale: d3Scale.scaleSequential,
    scheme: invertColorScale(d3ScaleChromatic.interpolateGreys),
    range: [0, 2250],
    empty: '#ccc',
  },
  terrestrialSlope: {
    scale: d3Scale.scaleSequential,
    scheme: d3ScaleChromatic.interpolateReds,
    range: [0, 90],
    empty: '#ccc',
  },
  droughtRiskWaterStress: {
    scale: d3Scale.scaleSequential,
    scheme: d3ScaleChromatic.interpolateReds,
    range: [0, 15],
    empty: '#ccc',
  },
  droughtRiskEpd: {
    scale: d3Scale.scaleSequential,
    scheme: d3ScaleChromatic.interpolateReds,
    range: [0, 2_500_000],
    empty: '#ccc',
  },
  droughtRiskEael: {
    scale: d3Scale.scaleSequential,
    scheme: invertColorScale(d3ScaleChromatic.interpolateInferno),
    range: [0, 1e9],
    empty: '#ccc',
  },
  droughtOptionsCost: {
    scale: d3Scale.scaleSequential,
    scheme: discardSides(d3ScaleChromatic.interpolateGreens, 0.2, 0.2),
    range: [0, 1_000_000_000],
    empty: '#ccc',
  },
  droughtOptionsPopulationProtected: {
    scale: d3Scale.scaleSequential,
    scheme: discardSides(d3ScaleChromatic.interpolateBlues, 0.2, 0.2),
    range: [0, 100_000],
    empty: '#ccc',
  },
  droughtOptionsNPVBenefit: {
    scale: d3Scale.scaleSequential,
    scheme: discardSides(d3ScaleChromatic.interpolateBlues, 0.2, 0.2),
    range: [0, 1_000_000_000],
    empty: '#ccc',
  },
  droughtOptionsBenefitCost: {
    scale: d3Scale.scaleSequential,
    scheme: invertColorScale(d3ScaleChromatic.interpolateViridis),
    range: [1, 1.5],
    empty: '#ccc',
  },
});
