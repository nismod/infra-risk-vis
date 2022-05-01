import { DataFilterExtension } from '@deck.gl/extensions';
import { TERRESTRIAL_LANDUSE_COLORS } from 'config/solutions/colors';
import { ViewLayer, FieldSpec } from 'lib/data-map/view-layers';
import { mvtLayer } from 'lib/deck/layers/base';
import { selector } from 'recoil';
import { sectionStyleValueState, sectionVisibilityState } from 'state/sections';
import { terrestrialFiltersState } from 'state/solutions/terrestrial-filters';
import { colorMap } from 'lib/color-map';
import { VECTOR_COLOR_MAPS } from 'config/color-maps';
import { featureProperty } from 'lib/deck/props/data-source';
import { dataColorMap } from 'lib/deck/props/color-map';
import { fillColor } from 'lib/deck/props/style';
import { Accessor } from 'lib/deck/props/getters';
import { LandUseOption, TerrestrialLocationFilterType } from 'config/solutions/domains';
import { truthyKeys } from 'lib/helpers';

function landuseColorMap(x: string) {
  return TERRESTRIAL_LANDUSE_COLORS[x].css;
}

export const terrestrialColorFnState = selector<Accessor<string>>({
  key: 'terrestrialColorFnState',
  get: ({ get }) => {
    const style = get(sectionStyleValueState('terrestrial'));

    if (style === 'landuse') {
      return landuseColorMap;
    } else if (style === 'elevation') {
      return colorMap(VECTOR_COLOR_MAPS.terrestrialElevation);
    } else if (style === 'slope') {
      return colorMap(VECTOR_COLOR_MAPS.terrestrialSlope);
    }
  },
});

export const terrestrialFieldSpecState = selector<FieldSpec>({
  key: 'terrestrialFieldSpecState',
  get: ({ get }) => {
    const style = get(sectionStyleValueState('terrestrial'));

    let field: string;

    if (style === 'landuse') {
      field = 'landuse_desc';
    } else if (style === 'elevation') {
      field = 'elevation_m';
    } else if (style === 'slope') {
      field = 'slope_degrees';
    }

    return {
      fieldGroup: 'properties',
      field,
    };
  },
});

function landuseFilterValue(p, landuseFilters: LandUseOption[]) {
  return landuseFilters.some((key) => p.landuse_desc === key) ? 1 : 0;
}

function locationFilterValue(p, locationFiltersKeys: TerrestrialLocationFilterType[]) {
  return locationFiltersKeys.every((key) => p[key] == true) ? 1 : 0;
}

export const terrestrialLayerState = selector<ViewLayer>({
  key: 'terrestrialLayerState',
  get: ({ get }) => {
    const showTerrestrial = get(sectionVisibilityState('terrestrial'));

    if (!showTerrestrial) {
      return null;
    }

    const filters = get(terrestrialFiltersState);
    const fieldSpec = get(terrestrialFieldSpecState);

    const dataFn = featureProperty(fieldSpec.field);
    const colorFn = get(terrestrialColorFnState);

    if (!colorFn || !dataFn) {
      return null;
    }

    const landuseFilterKeys = truthyKeys(filters.landuse_desc);
    const locationFilterKeys = truthyKeys(filters.location_filters);

    return {
      id: 'terrestrial',
      group: null,
      // interactionGroup: 'assets',
      fn: ({ deckProps }) => {
        return mvtLayer(
          deckProps,
          {
            data: '/vector/data/natural_terrestrial_combined.json',
            minZoom: 14,
            binary: false,
            filled: true,
            stroked: true,

            getLineColor: [250, 250, 250],
            getLineWidth: 2,
            lineWidthUnit: 'meters',
            getFilterValue: ({ properties }) => [
              landuseFilterValue(properties, landuseFilterKeys),
              properties.slope_degrees,
              properties.elevation_m,
              locationFilterValue(properties, locationFilterKeys),
            ],
            filterRange: [[1, 1], [...filters.slope_degrees], [...filters.elevation_m], [1, 1]],

            updateTriggers: {
              getFilterValue: [filters],
            },

            extensions: [new DataFilterExtension({ filterSize: 4 })],
          },
          fillColor(dataColorMap(dataFn, colorFn)),
        );
      },
    };
  },
});
