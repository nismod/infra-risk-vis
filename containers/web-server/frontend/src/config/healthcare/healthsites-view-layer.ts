import { fillColor, pointRadius } from '@/lib/deck/props/style';

import { assetViewLayer } from '../assets/asset-view-layer';
import { assetDataAccessFunction } from '../assets/data-access';
import { COLORS } from '../colors';

export function healthsitesViewLayer() {
  return assetViewLayer(
    'healthsites',
    {
      group: 'healthsites',
      spatialType: 'vector',
      interactionGroup: 'assets',
    },
    -1000,
    ({ zoom }) => [pointRadius(zoom), fillColor(COLORS.healthsites.deck)],
    assetDataAccessFunction('healthsites'),
  );
}
