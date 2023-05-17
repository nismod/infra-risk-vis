import { GeoJsonLayerProps } from 'deck.gl/typed';

import { Getter } from './getters';

export type ScaleLevel = 0 | 1 | 2 | 3;

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
  3: {
    getLineWidth: 10,
    lineWidthMinPixels: 1,
    lineWidthMaxPixels: 3,
  },
};

export function lineStyle(zoom, level: ScaleLevel = 2) {
  return {
    lineJointRounded: true,
    lineCapRounded: true,
    lineWidthUnit: 'meters',

    ...lineSizeLevels[level],
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
  3: { getPointRadius: 100, pointRadiusMinPixels: 0, pointRadiusMaxPixels: 3 },
};

export function pointRadius(zoom, level: ScaleLevel = 2): Partial<GeoJsonLayerProps> {
  return {
    pointRadiusUnits: 'meters',
    ...pointSizeLevels[level],
  };
}

export function iconSize(zoom: number, level: ScaleLevel = 2): Partial<GeoJsonLayerProps> {
  const { getPointRadius, pointRadiusMinPixels, pointRadiusMaxPixels } = pointSizeLevels[level];

  return {
    iconSizeUnits: 'meters',
    getIconSize: getPointRadius * 2,
    iconSizeMinPixels: pointRadiusMinPixels * 2,
    iconSizeMaxPixels: pointRadiusMaxPixels * 2,
  };
}

export type Color = [number, number, number, number?];
export type GetColor = Getter<Color>;

export function vectorColor(type: 'fill' | 'stroke' | 'icon', getColor: GetColor) {
  let propName: string;
  if (type === 'fill') propName = 'getFillColor';
  else if (type === 'stroke') propName = 'getLineColor';
  else if (type === 'icon') propName = 'getIconColor'; // this will only work for GeoJsonLayer

  return {
    [propName]: getColor,
    updateTriggers: {
      [propName]: (getColor as any)?.updateTriggers ?? (typeof getColor === 'function' ? [] : undefined),
    },
  };
}

export const iconColor = (getColor: GetColor) => vectorColor('icon', getColor);
export const fillColor = (getColor: GetColor) => vectorColor('fill', getColor);
export const strokeColor = (getColor: GetColor) => vectorColor('stroke', getColor);

export function border(color = [255, 255, 255], minPixels = 1) {
  return {
    stroked: true,
    getLineColor: color,
    lineWidthMinPixels: minPixels,
  };
}
