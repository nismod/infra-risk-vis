import React from 'react';

import { fillColor, pointRadius } from '@/lib/deck/props/style';
import { makeColor } from '@/lib/helpers';

import { VectorHoverDescription } from '@/map/tooltip/VectorHoverDescription';

import { assetViewLayer } from '../assets/asset-view-layer';
import { assetDataAccessFunction } from '../assets/data-access';
import { AssetMetadata } from '../assets/metadata';

export const HEALTHSITES_COLOR = makeColor('#72dfda');

export const HEALTHSITES_METADATA: AssetMetadata = {
  type: 'circle',
  label: 'Health Sites',
  color: HEALTHSITES_COLOR.css,
};

export function healthsitesViewLayer() {
  return assetViewLayer({
    assetId: 'healthsites',
    metadata: {
      spatialType: 'vector',
      interactionGroup: 'assets',
    },
    customFn: ({ zoom }) => [pointRadius(zoom), fillColor(HEALTHSITES_COLOR.deck)],
    customDataAccessFn: assetDataAccessFunction('healthsites'),
    renderTooltip: (hover) => {
      const { label, color } = HEALTHSITES_METADATA;
      return React.createElement(VectorHoverDescription, {
        hoveredObject: hover,
        label,
        color,
        idValue: hover.target.feature.properties.osm_id,
      });
    },
  });
}
