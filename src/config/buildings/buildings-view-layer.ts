import { ViewLayer } from 'lib/data-map/view-layers';
import { COLORS } from 'config/colors';
import { assetViewLayer } from 'config/assets/asset-view-layer';
import { border, fillColor } from 'lib/deck/props/style';
import { assetDataAccessFunction } from 'config/assets/data-accessor';

enum BuildingType {
  commercial = 'commercial',
  industrial = 'industrial',
  institutional = 'institutional',
  mixed_use = 'mixed_use',
  other = 'other',
  recreation = 'recreation',
  residential = 'residential',
  resort = 'resort',
}

const buildingTypeLookup = {
  Commercial: BuildingType.commercial,
  Industrial: BuildingType.industrial,
  Institutional: BuildingType.institutional,
  'Mixed Use': BuildingType.mixed_use,
  Other: BuildingType.other,
  Recreation: BuildingType.recreation,
  Residential: BuildingType.residential,
  Resort: BuildingType.resort,
};

const buildingColor = {
  [BuildingType.commercial]: COLORS.buildings_commercial.deck,
  [BuildingType.industrial]: COLORS.buildings_industrial.deck,
  [BuildingType.institutional]: COLORS.buildings_institutional.deck,
  [BuildingType.mixed_use]: COLORS.buildings_mixed_use.deck,
  [BuildingType.other]: COLORS.buildings_other.deck,
  [BuildingType.recreation]: COLORS.buildings_recreation.deck,
  [BuildingType.residential]: COLORS.buildings_residential.deck,
  [BuildingType.resort]: COLORS.buildings_resort.deck,
};

export function buildingsViewLayer(): ViewLayer {
  return assetViewLayer(
    'buildings',
    {
      group: 'assets',
      spatialType: 'vector',
      interactionGroup: 'assets',
    },
    -1000,
    ({ zoom, styleParams }) => [
      { minZoom: 12 },
      border(COLORS.buildings_unknown.deck),
      fillColor((x) => buildingColor[buildingTypeLookup[x.properties.building_type]]),
    ],
    assetDataAccessFunction('buildings'),
  );
}
