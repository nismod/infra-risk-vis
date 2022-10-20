import GL from '@luma.gl/constants';
import React from 'react';

import { ViewLayer } from '@/lib/data-map/view-layers';
import { rasterTileLayer } from '@/lib/deck/layers/raster-tile-layer';
import { numFormat } from '@/lib/helpers';

import { RasterLegend } from '@/map/legend/RasterLegend';

import { SOURCES } from '../sources';

export const JRC_POPULATION_COLOR_MAP: { scheme: string; range: [number, number] } = {
  scheme: 'purd',
  range: [0, 1e4],
};

function getPopulationUrl() {
  return SOURCES.raster.getUrl({
    path: 'population',
    ...JRC_POPULATION_COLOR_MAP,
  });
}

export function jrcPopulationViewLayer(): ViewLayer {
  return {
    id: 'population',
    group: 'population',
    interactionGroup: 'population',
    spatialType: 'raster',
    fn({ deckProps, zoom }) {
      return rasterTileLayer(
        {
          textureParameters: {
            [GL.TEXTURE_MAG_FILTER]: zoom >= 7 ? GL.NEAREST : GL.LINEAR,
          },
        },
        deckProps,
        {
          data: getPopulationUrl(),
          refinementStrategy: 'best-available',
        },
      );
    },
    renderLegend() {
      return React.createElement(RasterLegend, {
        key: 'population',
        label: 'Population',
        colorMap: JRC_POPULATION_COLOR_MAP,
        getValueLabel: (x) => numFormat(x),
      });
    },
  };
}
