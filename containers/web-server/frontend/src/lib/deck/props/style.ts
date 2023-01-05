import { GeoJsonLayerProps } from 'deck.gl/typed';

import { Getter } from './getters';

type ScaleLevel = 0 | 1 | 2;

const lineSizeLevels: Record<
  ScaleLevel,
  {
    getLineWidth: number;
    lineWidthMinPixels: number;
    lineWidthMaxPixels: number;
  }
> = {
  0: {
    getLineWidth: 90,
    lineWidthMinPixels: 1,
    lineWidthMaxPixels: 7,
  },
  1: {
    getLineWidth: 60,
    lineWidthMinPixels: 1,
    lineWidthMaxPixels: 6,
  },
  2: {
    getLineWidth: 20,
    lineWidthMinPixels: 1,
    lineWidthMaxPixels: 5,
  },
};

export function lineStyle(zoom, level: ScaleLevel = 2) {
  return {
    lineJointRounded: true,
    lineCapRounded: true,
    lineWidthUnit: 'meters',

    ...lineSizeLevels[level],

    // widthScale: 2 ** (15 - zoom),
  };
}

const pointSizeLevels: Record<
  ScaleLevel,
  {
    getPointRadius: number;
    pointRadiusMinPixels: number;
    pointRadiusMaxPixels: number;
  }
> = {
  0: { getPointRadius: 1500, pointRadiusMinPixels: 3, pointRadiusMaxPixels: 6 },
  1: { getPointRadius: 1200, pointRadiusMinPixels: 2, pointRadiusMaxPixels: 4 },
  2: { getPointRadius: 300, pointRadiusMinPixels: 2, pointRadiusMaxPixels: 4 },
};

export function pointRadius(zoom, level: ScaleLevel = 2): Partial<GeoJsonLayerProps> {
  return {
    pointRadiusUnits: 'meters',
    ...pointSizeLevels[level],
  };
}

export type Color = [number, number, number, number?];
export type GetColor = Getter<Color>;

/**
 * Returns color with alpha set to new value.
 * Doesn't mutate the input color.
 * @param color (r,g,b,[a]) color to modify
 * @param alpha new alpha value
 */
export function setAlpha(color: Color, alpha: number): Color {
  const [r, g, b] = color;
  return [r, g, b, alpha];
}

export function vectorColor(type: 'fill' | 'stroke', getColor: GetColor) {
  let propName: string;
  if (type === 'fill') propName = 'getFillColor';
  else if (type === 'stroke') propName = 'getLineColor';

  return {
    [propName]: getColor,
    updateTriggers: {
      [propName]:
        (getColor as any)?.updateTriggers ?? (typeof getColor === 'function' ? [] : undefined),
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
