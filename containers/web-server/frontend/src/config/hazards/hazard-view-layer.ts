import GL from '@luma.gl/constants';

import { ViewLayer } from '@/lib/data-map/view-layers';
import { rasterTileLayer } from '@/lib/deck/layers/raster-tile-layer';

import { RASTER_COLOR_MAPS } from '@/config/color-maps';
import { HazardParams } from '@/config/hazards/domains';

import { HAZARD_SOURCE } from './source';

export function getHazardId<
  F extends string, //'fluvial' | 'surface' | 'coastal' | 'cyclone',
  RP extends number,
  RCP extends string,
  E extends number | string,
  C extends number | string,
>({ hazardType, returnPeriod, rcp, epoch, gcm }: { hazardType: F; returnPeriod: RP; rcp: RCP; epoch: E; gcm: C }) {
  return `${hazardType}__rp_${returnPeriod}__rcp_${rcp}__epoch_${epoch}__conf_${gcm}` as const;
}

export function hazardViewLayer(hazardType: string, hazardParams: HazardParams): ViewLayer {
  const magFilter = hazardType === 'cyclone' ? GL.NEAREST : GL.LINEAR;

  const { returnPeriod, rcp, epoch, gcm } = hazardParams;

  const deckId = getHazardId({ hazardType, returnPeriod, rcp, epoch, gcm });

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
          data: HAZARD_SOURCE.getDataUrl({ hazardType, hazardParams }, { scheme, range }),
          refinementStrategy: 'no-overlap',
        },
      );
    },
  };
}
