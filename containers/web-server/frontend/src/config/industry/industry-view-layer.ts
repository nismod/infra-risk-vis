import React from 'react';

import { InteractionTarget, VectorTarget } from '@/lib/data-map/interactions/use-interactions';
import { fillColor, pointRadius } from '@/lib/deck/props/style';
import { makeColorConfig, makeConfig } from '@/lib/helpers';

import { VectorHoverDescription } from '@/map/tooltip/VectorHoverDescription';
import { IndustryType } from '@/state/data-selection/industry';

import { assetViewLayer } from '../assets/asset-view-layer';
import { assetDataAccessFunction } from '../assets/data-access';
import { AssetMetadata } from '../assets/metadata';

export const INDUSTRY_COLORS = makeColorConfig({
  cement: '#e4cda9',
  steel: '#5b8cc3',
});

export const INDUSTRY_METADATA = makeConfig<AssetMetadata, IndustryType>([
  {
    id: 'cement',
    type: 'circle',
    label: 'Industry (Cement)',
    color: INDUSTRY_COLORS.cement.css,
  },
  {
    id: 'steel',
    type: 'circle',
    label: 'Industry (Steel)',
    color: INDUSTRY_COLORS.steel.css,
  },
]);

export function industryViewLayer(industry_type_id: IndustryType) {
  return assetViewLayer({
    assetId: industry_type_id,
    metadata: {
      spatialType: 'vector',
      interactionGroup: 'assets',
    },
    customFn: ({ zoom }) => [pointRadius(zoom), fillColor(INDUSTRY_COLORS[industry_type_id].deck)],
    customDataAccessFn: assetDataAccessFunction(industry_type_id),
    renderTooltip: (hover: InteractionTarget<VectorTarget>) => {
      const { label, color } = INDUSTRY_METADATA[industry_type_id];

      return React.createElement(VectorHoverDescription, {
        hoveredObject: hover,
        label,
        color,
        idValue: hover.target.feature.properties.uid,
      });
    },
  });
}
