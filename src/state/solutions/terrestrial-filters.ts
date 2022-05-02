import { LandUseOption, LAND_USE_VALUES, TerrestrialLocationFilterType } from 'config/solutions/domains';
import _ from 'lodash';
import { atom } from 'recoil';

export type TerrestrialLocationFilters = Record<TerrestrialLocationFilterType, boolean>;

export interface TerrestrialFilters {
  landuse_desc: Record<LandUseOption, boolean>;
  slope_degrees: [number, number];
  elevation_m: [number, number];
  location_filters: TerrestrialLocationFilters;
}

export const terrestrialFiltersState = atom<TerrestrialFilters>({
  key: 'terrestrialFiltersState',
  default: {
    landuse_desc: _.fromPairs<boolean>(LAND_USE_VALUES.map((v) => [v, true])) as Record<LandUseOption, boolean>,
    slope_degrees: [0, 90],
    elevation_m: [0, 2250],
    location_filters: {
      within_forest_100m: false,
      is_protected: false,
      is_proposed_protected: false,
      within_major_river_50m: false,
      within_large_stream_50m: false,
      within_headwater_stream_50m: false,
    },
  },
});
