import React from 'react';

import { InteractionTarget, VectorTarget } from '@/lib/data-map/interactions/use-interactions';
import { border, fillColor, pointRadius } from '@/lib/deck/props/style';
import { makeColorConfig, makeConfig } from '@/lib/helpers';

import { SimpleAssetDetails } from '@/details/features/asset-details';
import { VectorHoverDescription } from '@/map/tooltip/VectorHoverDescription';
import { IndustryType } from '@/state/data-selection/industry';

import { assetViewLayer } from '../assets/asset-view-layer';
import { assetDataAccessFunction } from '../assets/data-access';
import { AssetMetadata } from '../assets/metadata';
import { IndustryDetails } from './details';

export const INDUSTRY_COLORS = makeColorConfig({
  cement: '#e4cda9',
  steel: '#5b8cc3',
});

export const INDUSTRY_METADATA = makeConfig<AssetMetadata & { shortLabel: string }, IndustryType>([
  {
    id: 'cement',
    type: 'circle',
    label: 'Industry (Cement)',
    shortLabel: 'Cement',
    color: INDUSTRY_COLORS.cement.css,
  },
  {
    id: 'steel',
    type: 'circle',
    label: 'Industry (Steel)',
    shortLabel: 'Steel',
    color: INDUSTRY_COLORS.steel.css,
  },
]);

export function industryViewLayer(industry_type_id: IndustryType) {
  const { label, color } = INDUSTRY_METADATA[industry_type_id];

  return assetViewLayer({
    assetId: industry_type_id,
    metadata: {
      spatialType: 'vector',
      interactionGroup: 'assets',
    },
    customFn: ({ zoom }) => [
      pointRadius(zoom, 1),
      fillColor(INDUSTRY_COLORS[industry_type_id].deck),
      border([100, 100, 100]),
    ],
    customDataAccessFn: assetDataAccessFunction(industry_type_id),
    renderTooltip: (hover: InteractionTarget<VectorTarget>) => {
      return React.createElement(VectorHoverDescription, {
        hoveredObject: hover,
        label,
        color,
        idValue: hover.target.feature.properties.uid,
      });
    },
    renderDetails(selection: InteractionTarget<VectorTarget>) {
      const feature = selection.target.feature;

      return React.createElement(SimpleAssetDetails, {
        feature,
        label,
        color,
        detailsComponent: IndustryDetails,
      });
    },
  });
}
