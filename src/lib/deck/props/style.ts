import { Getter } from './getters';

export function lineStyle(zoom) {
  return {
    getLineWidth: 15,
    lineWidthUnit: 'meters',
    lineWidthMinPixels: 1,
    lineWidthMaxPixels: 5,
    lineJointRounded: true,
    lineCapRounded: true,

    // widthScale: 2 ** (15 - zoom),
  };
}

export function pointRadius(zoom) {
  return {
    getPointRadius: 20,
    pointRadiusUnit: 'meters',
    pointRadiusMinPixels: 3,
    pointRadiusMaxPixels: 10,
    // radiusScale: 2 ** (15 - zoom),
  };
}

export type Color = [number, number, number, number?];
export type GetColor = Getter<Color>;

export function vectorColor(type: 'fill' | 'stroke', getColor: GetColor) {
  let propName: string;
  if (type === 'fill') propName = 'getFillColor';
  else if (type === 'stroke') propName = 'getLineColor';

  return {
    [propName]: getColor,
    updateTriggers: {
      [propName]: (getColor as any)?.updateTriggers ?? (typeof getColor === 'function' ? [] : undefined),
    },
  };
}

export const fillColor = (getColor: GetColor) => vectorColor('fill', getColor);
export const strokeColor = (getColor: GetColor) => vectorColor('stroke', getColor);

export function border(color = [255, 255, 255]) {
  return {
    stroked: true,
    getLineColor: color,
    lineWidthMinPixels: 1,
  };
}
