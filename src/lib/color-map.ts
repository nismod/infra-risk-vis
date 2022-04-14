import * as d3 from 'd3-scale';

export function colorMap(scale: (t: number, n?: number) => string, range: number[], empty: string) {
  const scaleFn = d3.scaleSequential<string>(range, scale);

  return (value) => (value == null ? empty : scaleFn(value));
}
