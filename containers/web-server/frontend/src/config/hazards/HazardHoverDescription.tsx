import { FC } from 'react';

import { InteractionTarget, RasterTarget } from '@/lib/data-map/interactions/use-interactions';

import { HAZARDS_METADATA, HAZARD_COLOR_MAPS } from '@/config/hazards/metadata';
import { RasterHoverDescription } from '@/map/tooltip/RasterHoverDescription';

export const HazardHoverDescription: FC<{ hoveredObject: InteractionTarget<RasterTarget> }> = ({ hoveredObject }) => {
  const { color } = hoveredObject.target;

  const {
    viewLayer: {
      id,
      params: { hazardType },
    },
  } = hoveredObject;
  const { label, dataUnit, fractionDigits } = HAZARDS_METADATA[id];
  const { scheme, range } = HAZARD_COLOR_MAPS[hazardType];

  return (
    <RasterHoverDescription
      color={color}
      scheme={scheme}
      range={range}
      label={`${label}`}
      formatValue={(x) =>
        x != null
          ? `${x.toLocaleString(undefined, {
              maximumFractionDigits: fractionDigits,
            })}${dataUnit}`
          : ''
      }
    />
  );
};
