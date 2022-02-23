import GL from '@luma.gl/constants';
import { HazardParams } from 'config/hazards/domains';

import { rasterTileLayer } from 'lib/deck-layers/raster-tile-layer';
import { ViewLayer } from 'lib/data-map/view-layers';

import { RASTER_COLOR_MAPS } from '../color-maps';

export function getHazardId<
  F extends string, //'fluvial' | 'surface' | 'coastal' | 'cyclone',
  RP extends number,
  RCP extends string,
  E extends number,
  C extends number | string,
>({
  hazardType,
  returnPeriod,
  rcp,
  epoch,
  confidence,
}: {
  hazardType: F;
  returnPeriod: RP;
  rcp: RCP;
  epoch: E;
  confidence: C;
}) {
  return `${hazardType}__rp_${returnPeriod}__rcp_${rcp}__epoch_${epoch}__conf_${confidence}` as const;
}

export function hazardViewLayer(hazardType: string, hazardParams: HazardParams): ViewLayer {
  const magFilter = hazardType === 'cyclone' ? GL.NEAREST : GL.LINEAR;

  const { returnPeriod, rcp, epoch, confidence } = hazardParams;
  const sanitisedRcp = rcp.replace('.', 'x');

  const deckId = getHazardId({ hazardType, returnPeriod, rcp, epoch, confidence }); //`hazard_${hazardType}_${returnPeriod}`;

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
          opacity: hazardType === 'cyclone' ? 0.6 : 1,
        },
        deckProps,
        {
          id: `${hazardType}@${deckId}`, // follow the convention viewLayerId@deckLayerId
          data: `/raster/singleband/${hazardType}/${returnPeriod}/${sanitisedRcp}/${epoch}/${confidence}/{z}/{x}/{y}.png?colormap=${scheme}&stretch_range=[${range[0]},${range[1]}]`,
          refinementStrategy: 'no-overlap',
        },
      );
    },
  };
}
