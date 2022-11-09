import GL from '@luma.gl/constants';
import React from 'react';

import { InteractionTarget, RasterTarget } from '@/lib/data-map/interactions/use-interactions';
import { ViewLayer } from '@/lib/data-map/view-layers';
import { rasterTileLayer } from '@/lib/deck/layers/raster-tile-layer';

import { HazardLegend } from '@/map/legend/content/HazardLegend';

import { HazardHoverDescription } from './HazardHoverDescription';
import { HAZARD_COLOR_MAPS } from './metadata';
import { getHazardDataPath, getHazardDataUrl } from './source';

export function getHazardId(hazardType: string, hazardParams: any) {
  return getHazardDataPath({ hazardType, metric: 'occurrence', hazardParams });
}

export function hazardViewLayer(hazardType: string, hazardParams: any): ViewLayer {
  const magFilter = hazardType === 'cyclone' ? GL.NEAREST : GL.LINEAR;

  const id = hazardType;
  const deckId = getHazardId(hazardType, hazardParams);

  return {
    id,
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
          opacity: hazardType === 'cyclone' ? 0.6 : 1,

          // TODO: tweak transparentColor to tweak border color / transparent layer tint
          /*
          transparentColor: [128, 128, 128, 0],
          */
        },
        deckProps,
        {
          id: `${id}@${deckId}`, // follow the convention viewLayerId@deckLayerId
          data: getHazardDataUrl({ hazardType, metric: 'occurrence', hazardParams }, { scheme, range }),
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
    renderTooltip(hover: InteractionTarget<RasterTarget>) {
      return React.createElement(HazardHoverDescription, {
        hoveredObject: hover,
      });
    },
  };
}
