import GL from '@luma.gl/constants';
import React from 'react';

import { ViewLayer } from '@/lib/data-map/view-layers';
import { rasterTileLayer } from '@/lib/deck/layers/raster-tile-layer';

import { RASTER_COLOR_MAPS } from '@/config/color-maps';
import { HazardLegend } from '@/map/legend/content/HazardLegend';

import { HAZARD_SOURCE } from './source';

export function getHazardId({ hazardType, hazardParams }: { hazardType: string; hazardParams: any }) {
  if (hazardType === 'earthquake') {
    const { returnPeriod, medium } = hazardParams;

    return `${hazardType}__rp_${returnPeriod}__medium_${medium}`;
  } else {
    const { returnPeriod, rcp, epoch, gcm } = hazardParams;

    return `${hazardType}__rp_${returnPeriod}__rcp_${rcp}__epoch_${epoch}__conf_${gcm}`;
  }
}

export function hazardViewLayer(hazardType: string, hazardParams: any): ViewLayer {
  const magFilter = hazardType === 'cyclone' ? GL.NEAREST : GL.LINEAR;

  const deckId = getHazardId({ hazardType, hazardParams });

  return {
    id: hazardType,
    group: 'hazards',
    spatialType: 'raster',
    interactionGroup: 'hazards',
    params: { hazardType, hazardParams },
    fn: ({ deckProps, zoom }) => {
      const { scheme, range } = RASTER_COLOR_MAPS[hazardType];

      return rasterTileLayer(
        {
          textureParameters: {
            [GL.TEXTURE_MAG_FILTER]: magFilter,
            // [GL.TEXTURE_MAG_FILTER]: zoom < 12 ? GL.NEAREST : GL.NEAREST_MIPMAP_LINEAR,
          },
          opacity: hazardType === 'cyclone' || hazardType === 'earthquake' ? 0.6 : 1,
        },
        deckProps,
        {
          id: `${hazardType}@${deckId}`, // follow the convention viewLayerId@deckLayerId
          data: HAZARD_SOURCE.getDataUrl({ hazardType, hazardParams }, { scheme, range }),
          refinementStrategy: 'no-overlap',
        },
      );
    },
    renderLegend() {
      return React.createElement(HazardLegend, {
        key: hazardType,
        viewLayer: this,
      });
    },
  };
}
