import GL from '@luma.gl/constants';
import _ from 'lodash';
import React from 'react';

import { InteractionTarget, RasterTarget } from '@/lib/data-map/interactions/use-interactions';
import { ViewLayer } from '@/lib/data-map/view-layers';
import { rasterTileLayer } from '@/lib/deck/layers/raster-tile-layer';

import { RasterBaseHover, serializeColor } from '@/map/tooltip/RasterBaseHover';

import landCoverLegend from './land-cover-legend.json';

const landCoverColorMap = _.map(landCoverLegend, ({ color: rgba }, code) => ({
  value: parseInt(code, 10),
  color: serializeColor(rgba as any),
}));

const landCoverLabels = Object.fromEntries(_.map(landCoverLegend, ({ name }, code) => [parseInt(code, 10), name]));

export function landCoverViewLayer(): ViewLayer {
  return {
    id: 'land_cover',
    spatialType: 'raster',
    interactionGroup: 'raster_assets',
    fn: ({ deckProps }) =>
      rasterTileLayer(
        {
          textureParameters: {
            [GL.TEXTURE_MAG_FILTER]: GL.LINEAR,
          },
        },
        deckProps,
        {
          data: '/api/tiles/land_cover/{z}/{x}/{y}.png?colormap=explicit',
          refinementStrategy: 'no-overlap',
        },
      ),
    renderTooltip(hover: InteractionTarget<RasterTarget>) {
      return React.createElement(RasterBaseHover, {
        colorMapValues: landCoverColorMap,
        color: hover.target.color,
        label: 'Land Cover',
        formatValue: (x) => landCoverLabels[x],
      });
    },
  };
}
