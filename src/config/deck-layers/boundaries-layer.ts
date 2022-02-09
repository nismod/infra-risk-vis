import { MVTLayer } from 'deck.gl';

import { border } from 'lib/deck-layers/utils';

import { BackgroundName } from '../backgrounds';

export type BoundaryLevel = 'parish' | 'community' | 'enumeration';

interface BoundaryConfig {
  fieldName: string;
  minZoom: number;
  showLabels: boolean;
}

export const boundaryConfig: Record<BoundaryLevel, BoundaryConfig> = {
  parish: {
    fieldName: 'PARISH',
    minZoom: 9,
    showLabels: true,
  },
  community: {
    fieldName: 'COMMUNITY',
    minZoom: 13,
    showLabels: false,
  },
  enumeration: {
    fieldName: 'ED',
    minZoom: 13,
    showLabels: false,
  },
};

export function boundariesLayer(level: BoundaryLevel) {
  return new MVTLayer(
    {
      id: `boundaries_${level}`,
      data: `/vector/data/boundaries_${level}.json`,
      loadOptions: {
        mvt: {
          layers: ['polygons'],
        },
      },
      binary: true,
      filled: true,
      getFillColor: [255, 255, 255, 0],
      pickable: true,
      stroked: true,
      refinementStrategy: 'best-available',
      lineWidthUnits: 'pixels',
    } as any,
    border([150, 150, 150, 255]) as any,
  );
}

export function boundaryLabelsLayer(level: BoundaryLevel, background: BackgroundName) {
  const config = boundaryConfig[level];

  const color = background === 'satellite' ? [240, 240, 240, 255] : [90, 90, 90, 255];
  return (
    config.showLabels &&
    new MVTLayer({
      id: `boundaries_${level}-text`,
      data: `/vector/data/boundaries_${level}_labels.json`,
      loadOptions: {
        mvt: {
          layers: ['labels'],
        },
      },
      binary: false,
      minZoom: config.minZoom,
      pointType: 'text',
      getText: (f) => f.properties[config.fieldName],
      getTextSize: 24,
      getTextColor: color,
      textFontFamily: 'Arial',
      textFontWeight: 'bold',
      getPolygonOffset: ({ layerIndex }) => [0, -layerIndex * 100 - 2000],

      // won't work before deck.gl v8.7.0 is released (textFontSettings isn't mapped correctly)
      // see https://github.com/visgl/deck.gl/pull/6336
      //
      // textOutlineColor: [255, 255, 255, 255],
      // textOutlineWidth: 1,
      // textFontSettings: {
      //   sdf: true,
      // },
    } as any)
  );
}
