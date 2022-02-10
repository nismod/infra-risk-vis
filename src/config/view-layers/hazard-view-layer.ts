import GL from '@luma.gl/constants';

import { rasterTileLayer } from 'lib/deck-layers/raster-tile-layer';
import { ViewLayer } from 'lib/view-layers';

import { RASTER_COLOR_MAPS } from '../color-maps';
import { getHazardId } from '../layers';

export function hazardViewLayer(hazardType, returnPeriod, rcp, epoch, confidence): ViewLayer {
  const id = getHazardId({ hazardType, returnPeriod, rcp, epoch, confidence }); //`hazard_${hazardType}_${returnPeriod}`;

  const magFilter = hazardType === 'cyclone' ? GL.NEAREST : GL.LINEAR;

  const sanitisedRcp = rcp.replace('.', 'x');

  return {
    id,
    spatialType: 'raster',
    interactionGroup: 'hazards',
    dataParams: { hazardType, returnPeriod, rcp, epoch, confidence },
    fn: ({ props, zoom, params: { hazardType, returnPeriod, rcp, epoch, confidence } }) => {
      const { scheme, range } = RASTER_COLOR_MAPS[hazardType];

      return rasterTileLayer(
        {
          textureParameters: {
            [GL.TEXTURE_MAG_FILTER]: magFilter,
            // [GL.TEXTURE_MAG_FILTER]: zoom < 12 ? GL.NEAREST : GL.NEAREST_MIPMAP_LINEAR,
          },
        },
        props,
        {
          id,
          data: `/raster/singleband/${hazardType}/${returnPeriod}/${sanitisedRcp}/${epoch}/${confidence}/{z}/{x}/{y}.png?colormap=${scheme}&stretch_range=[${range[0]},${range[1]}]`,
          refinementStrategy: 'no-overlap',
        },
      );
    },
  };
}
