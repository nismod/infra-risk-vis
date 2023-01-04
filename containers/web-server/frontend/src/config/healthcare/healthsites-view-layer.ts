import React from 'react';

import { border, fillColor, pointRadius } from '@/lib/deck/props/style';
import { makeColor } from '@/lib/helpers';

import { SimpleAssetDetails } from '@/details/features/asset-details';
import { VectorHoverDescription } from '@/map/tooltip/VectorHoverDescription';

import { assetViewLayer } from '../assets/asset-view-layer';
import { assetDataAccessFunction } from '../assets/data-access';
import { AssetMetadata } from '../assets/metadata';
import { HealthsiteDetails } from './details';

export const HEALTHSITES_COLOR = makeColor('#72dfda');

export const HEALTHSITES_METADATA: AssetMetadata = {
  type: 'circle',
  label: 'Healthcare',
  color: HEALTHSITES_COLOR.css,
};

export function healthsitesViewLayer() {
  const { label, color } = HEALTHSITES_METADATA;

  return assetViewLayer({
    assetId: 'healthsites',
    metadata: {
      spatialType: 'vector',
      interactionGroup: 'assets',
    },
    customFn: ({ zoom }) => [
      pointRadius(zoom),
      fillColor(HEALTHSITES_COLOR.deck),
      border([255, 255, 255]),
    ],
    customDataAccessFn: assetDataAccessFunction('healthsites'),
    renderTooltip: (hover) => {
      return React.createElement(VectorHoverDescription, {
        hoveredObject: hover,
        label,
        color,
        idValue: hover.target.feature.properties.osm_id,
      });
    },
    renderDetails(selection) {
      const feature = selection.target.feature;

      return React.createElement(SimpleAssetDetails, {
        label,
        color,
        feature,
        detailsComponent: HealthsiteDetails,
      });
    },
  });
}
