import React from 'react';

import { InteractionTarget, RasterTarget } from '@/lib/data-map/interactions/use-interactions';
import { ViewLayer } from '@/lib/data-map/view-layers';
import { rasterTileLayer } from '@/lib/deck/layers/raster-tile-layer';
import { Subset } from '@/lib/helpers';

import { RasterColorMap, RasterLegend } from '@/map/legend/RasterLegend';
import { RasterHoverDescription } from '@/map/tooltip/RasterHoverDescription';

import { HAZARDS_METADATA, HazardType } from '../metadata';
import { getHazardDataPath, getHazardDataUrl } from '../source';

export type ExposureSource = Subset<HazardType, 'extreme_heat' | 'drought'>;

export const EXPOSURE_COLOR_MAPS: Record<ExposureSource, RasterColorMap> = {
  extreme_heat: {
    scheme: 'reds',
    range: [0, 1000],
  },
  drought: {
    scheme: 'oranges',
    range: [0, 1000],
  },
};

function numFormatWhole(value: number) {
  return `${value.toLocaleString(undefined, { maximumFractionDigits: 0 })}`;
}

export function exposureViewLayer(hazardType: ExposureSource, hazardParams: any): ViewLayer {
  const id = `${hazardType}_exposure`;
  const deckId = getHazardDataPath({ hazardType, hazardParams, metric: 'exposure' });

  let { label: hazardLabel } = HAZARDS_METADATA[hazardType];
  const colorMap = EXPOSURE_COLOR_MAPS[hazardType];

  const label = `Pop. Exposed (${hazardLabel})`;

  return {
    id,
    spatialType: 'raster',
    interactionGroup: 'hazards',
    params: { hazardType, hazardParams },
    fn: ({ deckProps, zoom }) => {
      return rasterTileLayer({}, deckProps, {
        id: `${id}@${deckId}`, // follow the convention viewLayerId@deckLayerId
        data: getHazardDataUrl({ hazardType, metric: 'exposure', hazardParams }, colorMap),
        refinementStrategy: 'no-overlap',
      });
    },
    renderLegend() {
      return React.createElement(RasterLegend, {
        label,
        colorMap,
        getValueLabel: numFormatWhole,
      });
    },
    renderTooltip(hover: InteractionTarget<RasterTarget>) {
      return React.createElement(RasterHoverDescription, {
        ...colorMap,
        color: hover.target.color,
        label,
        formatValue: numFormatWhole,
      });
    },
  };
}
