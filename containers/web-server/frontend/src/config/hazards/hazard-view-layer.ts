import GL from '@luma.gl/constants';
import React from 'react';

import { ViewLayer } from '@/lib/data-map/view-layers';
import { rasterTileLayer } from '@/lib/deck/layers/raster-tile-layer';

import { HazardLegend } from '@/map/legend/content/HazardLegend';

import { HAZARD_COLOR_MAPS } from './metadata';
import { HAZARD_SOURCE } from './source';

export function getHazardId({ hazardType, hazardParams }: { hazardType: string; hazardParams: any }) {
  if (hazardType === 'earthquake') {
    const { rp, medium } = hazardParams;

    return `${hazardType}__rp_${rp}__medium_${medium}`;
  } else if (hazardType === 'drought') {
    const { rcp, epoch, gcm } = hazardParams;

    return `${hazardType}__rcp_${rcp}__epoch_${epoch}__gcm_${gcm}`;
  } else {
    const { rp, rcp, epoch, gcm } = hazardParams;

    return `${hazardType}__rp_${rp}__rcp_${rcp}__epoch_${epoch}__gcm_${gcm}`;
  }
}

export function hazardViewLayer(hazardType: string, hazardParams: any): ViewLayer {
  const magFilter = hazardType === 'cyclone' ? GL.NEAREST : GL.LINEAR;

  const deckId = getHazardId({ hazardType, hazardParams });

  return {
    id: hazardType,
    spatialType: 'raster',
    interactionGroup: 'hazards',
    params: { hazardType, hazardParams },
    fn: ({ deckProps, zoom }) => {
      const { scheme, range } = HAZARD_COLOR_MAPS[hazardType];

      return rasterTileLayer(
        {
          textureParameters: {
            [GL.TEXTURE_MAG_FILTER]: magFilter,
          },
          opacity: hazardType === 'cyclone' || hazardType === 'earthquake' ? 0.6 : 1,

          // TODO: tweak transparentColor to tweak border color / transparent layer tint
          /*
          transparentColor: [128, 128, 128, 0],
          */
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
