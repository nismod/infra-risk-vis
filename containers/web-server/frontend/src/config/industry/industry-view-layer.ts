import { fillColor, pointRadius } from '@/lib/deck/props/style';

import { assetViewLayer } from '../assets/asset-view-layer';
import { assetDataAccessFunction } from '../assets/data-access';
import { COLORS } from '../colors';

const industryColor = {
  cement: COLORS.industry_cement.deck,
  steel: COLORS.industry_steel.deck,
};

export function industryViewLayer(industry_type_id) {
  return assetViewLayer(
    industry_type_id,
    {
      group: 'industry',
      spatialType: 'vector',
      interactionGroup: 'assets',
    },
    -1000,
    ({ zoom }) => [pointRadius(zoom), fillColor(industryColor[industry_type_id])],
    assetDataAccessFunction(industry_type_id),
  );
}
