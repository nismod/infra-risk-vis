import React from 'react';

import { InteractionTarget, RasterTarget } from '@/lib/data-map/interactions/use-interactions';
import { ViewLayer } from '@/lib/data-map/view-layers';
import { rasterTileLayer } from '@/lib/deck/layers/raster-tile-layer';
import { makeValueFormat, nullFormat } from '@/lib/formats';
import { toLabelLookup } from '@/lib/helpers';

import { RasterLegend } from '@/map/legend/RasterLegend';
import { RasterHoverDescription } from '@/map/tooltip/RasterHoverDescription';

import { SOURCES } from '../sources';
import { BUILDING_DENSITY_COLORMAPS, BUILDING_DENSITY_TYPE_LABELS, BuildingDensityType } from './metadata';

const typeLabelLookup = toLabelLookup(BUILDING_DENSITY_TYPE_LABELS);

export function buildingDensityLayer(type: BuildingDensityType): ViewLayer {
  const colorMap = BUILDING_DENSITY_COLORMAPS[type];

  const label = `Built-up Surface (${typeLabelLookup[type]})`;

  const formatFn = makeValueFormat(
    (x) => (
      <>
        {x}m<sup>2</sup>
      </>
    ),
    { maximumFractionDigits: 0 },
  );

  const formatValue = nullFormat(formatFn);

  return {
    id: 'buildings',
    interactionGroup: 'raster_assets',
    spatialType: 'raster',
    fn: ({ deckProps }) => {
      return rasterTileLayer(
        {
          transparentColor: [255, 255, 255, 0],
        },
        deckProps,
        {
          data: SOURCES.raster.getUrl({
            path: `buildings/${type}`,
            ...colorMap,
          }),
        },
      );
    },
    renderLegend() {
      return React.createElement(RasterLegend, {
        label,
        colorMap,
        getValueLabel: formatValue,
      });
    },
    renderTooltip(hoveredObject: InteractionTarget<RasterTarget>) {
      return React.createElement(RasterHoverDescription, {
        ...colorMap,
        color: hoveredObject.target.color,
        label,
        formatValue,
      });
    },
  };
}
