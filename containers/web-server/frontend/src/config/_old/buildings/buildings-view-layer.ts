import { ViewLayer } from '@/lib/data-map/view-layers';
import { border, fillColor } from '@/lib/deck/props/style';

import { assetViewLayer } from '@/config/assets/asset-view-layer';
import { assetDataAccessFunction } from '@/config/assets/data-access';
import { BuildingType } from '@/state/data-selection/_old/buildings';

import { BUILDING_COLORS } from './metadata';

const buildingColor = {
  buildings_commercial: BUILDING_COLORS.commercial,
  buildings_industrial: BUILDING_COLORS.industrial,
  buildings_institutional: BUILDING_COLORS.institutional,
  buildings_mixed: BUILDING_COLORS.mixed,
  buildings_other: BUILDING_COLORS.other,
  buildings_recreation: BUILDING_COLORS.recreation,
  buildings_residential: BUILDING_COLORS.residential,
  buildings_resort: BUILDING_COLORS.resort,
};

export function buildingsViewLayer(building_type_id: BuildingType): ViewLayer {
  return assetViewLayer({
    assetId: building_type_id,
    metadata: {
      spatialType: 'vector',
      interactionGroup: 'assets',
    },
    customFn: () => [
      { minZoom: 12 },
      border(BUILDING_COLORS.unknown.deck),
      fillColor(buildingColor[building_type_id].deck),
    ],
    customDataAccessFn: assetDataAccessFunction(building_type_id),
  });
}
