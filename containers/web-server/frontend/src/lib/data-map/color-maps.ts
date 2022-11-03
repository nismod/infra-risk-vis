import * as d3Scale from 'd3-scale';
import * as d3ScaleChromatic from 'd3-scale-chromatic';

export { d3Scale, d3ScaleChromatic };

export function invertColorScale<T>(colorScale: (t: number) => T) {
  return (i: number) => colorScale(1 - i);
}

export function discardSides<T>(interpolator: (t: number) => T, cutStart: number, cutEnd: number = 0) {
  return (i: number) => {
    const t = i * (1 - cutStart - cutEnd) + cutStart;
    return interpolator(t);
  };
}
