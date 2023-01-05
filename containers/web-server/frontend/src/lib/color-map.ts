import * as d3Array from 'd3-array';

import { ColorSpec } from './data-map/view-layers';

export function colorScaleFn({ scale, range, scheme }: ColorSpec) {
  return scale(range, scheme);
}

export function colorScaleValues(colorSpec: ColorSpec, n: number) {
  const scaleFn = colorScaleFn(colorSpec);
  const [rangeMin, rangeMax] = colorSpec.range;
  return d3Array.ticks(rangeMin, rangeMax, n).map((x) => ({ value: x, color: scaleFn(x) }));
}

export function colorMap(colorSpec: ColorSpec) {
  const scaleFn = colorScaleFn(colorSpec);

  return (value) =>
    value == null || (colorSpec.zeroIsEmpty && value === 0) ? colorSpec.empty : scaleFn(value);
}
