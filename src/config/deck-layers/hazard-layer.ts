import GL from '@luma.gl/constants';

import { RASTER_COLOR_MAPS } from '../color-maps';
import { getHazardId } from '../layers';

import { rasterTileLayer } from './raster-tile-layer';

export function hazardDeckLayer(hazardType, returnPeriod, rcp, epoch, confidence) {
  const id = getHazardId({ hazardType, returnPeriod, rcp, epoch, confidence }); //`hazard_${hazardType}_${returnPeriod}`;

  const magFilter = hazardType === 'cyclone' ? GL.NEAREST : GL.LINEAR;
  const refinementStrategy = hazardType === 'cyclone' ? 'best-available' : 'no-overlap';

  const sanitisedRcp = rcp.replace('.', 'x');

  return {
    id,
    spatialType: 'raster',
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
          refinementStrategy,
        },
      );
    },
  };
}