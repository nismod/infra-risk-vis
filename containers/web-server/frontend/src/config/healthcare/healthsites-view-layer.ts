import { fillColor, pointRadius } from '@/lib/deck/props/style';

import { assetViewLayer } from '../assets/asset-view-layer';
import { assetDataAccessFunction } from '../assets/data-access';
import { COLORS } from '../colors';

export function healthsitesViewLayer() {
  return assetViewLayer({
    assetId: 'healthsites',
    metadata: {
      spatialType: 'vector',
      interactionGroup: 'assets',
    },
    customFn: ({ zoom }) => [pointRadius(zoom), fillColor(COLORS.healthsites.deck)],
    customDataAccessFn: assetDataAccessFunction('healthsites'),
  });
}
