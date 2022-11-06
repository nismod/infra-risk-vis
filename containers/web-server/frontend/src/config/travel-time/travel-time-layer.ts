import { capitalize } from 'lodash';
import React from 'react';

import { InteractionTarget, RasterTarget } from '@/lib/data-map/interactions/use-interactions';
import { rasterTileLayer } from '@/lib/deck/layers/raster-tile-layer';

import { RasterColorMap, RasterLegend } from '@/map/legend/RasterLegend';
import { RasterHoverDescription } from '@/map/tooltip/RasterHoverDescription';

import { SOURCES } from '../sources';

export const TRAVELTIME_TYPES = ['motorized', 'walking'] as const;

export type TraveltimeType = typeof TRAVELTIME_TYPES[number];

export const TRAVELTIME_VALUE_LABELS = TRAVELTIME_TYPES.map((x) => ({ value: x, label: capitalize(x) }));

const TRAVELTIME_COLORMAP: RasterColorMap = {
  //TODO
  scheme: 'rdbu_r',
  range: [0, 240],
};

export function travelTimeViewLayer(type: TraveltimeType) {
  const id = `traveltime_to_healthcare_${type}`;
  const label = `Travel Time to Healthcare (${type})`;

  const formatValue = (x: number) => (x == null ? '-' : `${x.toLocaleString(undefined, { maximumFractionDigits: 1 })}`);

  return {
    id,
    spatialType: 'raster',
    interactionGroup: 'raster_assets',
    fn({ deckProps }) {
      return rasterTileLayer({}, deckProps, {
        data: SOURCES.raster.getUrl({
          path: `traveltime_to_healthcare/${type}`,
          ...TRAVELTIME_COLORMAP,
        }),
        refinementStrategy: 'no-overlap',
      });
    },
    renderLegend() {
      return React.createElement(RasterLegend, {
        label,
        colorMap: TRAVELTIME_COLORMAP,
        getValueLabel: formatValue,
      });
    },
    renderTooltip(hover: InteractionTarget<RasterTarget>) {
      const { color } = hover.target;
      return React.createElement(RasterHoverDescription, {
        ...TRAVELTIME_COLORMAP,
        color,
        label,
        formatValue,
      });
    },
  };
}
