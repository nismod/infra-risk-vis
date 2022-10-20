import { FC } from 'react';

import { InteractionTarget, RasterTarget } from '@/lib/data-map/interactions/use-interactions';
import { numFormat } from '@/lib/helpers';

import { JRC_POPULATION_COLOR_MAP } from '@/config/population/population-view-layer';

import { RasterHoverDescription } from './RasterHoverDescription';

export const PopulationHoverDescription: FC<{ hoveredObject: InteractionTarget<RasterTarget> }> = ({
  hoveredObject,
}) => {
  const { color } = hoveredObject.target;
  return (
    <RasterHoverDescription
      color={color}
      {...JRC_POPULATION_COLOR_MAP}
      label="Population"
      formatValue={(x) => numFormat(x)}
    />
  );
};
