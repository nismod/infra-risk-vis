import { FC } from 'react';

import { InteractionTarget, RasterTarget } from '@/lib/data-map/interactions/use-interactions';

import { RASTER_COLOR_MAPS } from '@/config/color-maps';
import { HAZARDS_METADATA } from '@/config/hazards/metadata';

import { RasterHoverDescription } from './RasterHoverDescription';

export const HazardHoverDescription: FC<{ hoveredObject: InteractionTarget<RasterTarget> }> = ({ hoveredObject }) => {
  const { color } = hoveredObject.target;

  const {
    viewLayer: {
      id,
      params: { hazardType },
    },
  } = hoveredObject;
  const { label, dataUnit } = HAZARDS_METADATA[id];
  const { scheme, range } = RASTER_COLOR_MAPS[hazardType];

  return (
    <RasterHoverDescription
      color={color}
      scheme={scheme}
      range={range}
      label={`${label}`}
      formatValue={(x) => (x != null ? `${x.toFixed(1)}${dataUnit}` : '')}
    />
  );
};
