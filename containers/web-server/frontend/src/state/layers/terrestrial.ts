import { DataFilterExtension } from '@deck.gl/extensions';
import { selector } from 'recoil';

import { colorMap } from '@/lib/color-map';
import { ColorSpec, FieldSpec, ViewLayer } from '@/lib/data-map/view-layers';
import { selectableMvtLayer } from '@/lib/deck/layers/selectable-mvt-layer';
import { dataColorMap } from '@/lib/deck/props/color-map';
import { featureProperty } from '@/lib/deck/props/data-source';
import { Accessor } from '@/lib/deck/props/getters';
import { border, fillColor } from '@/lib/deck/props/style';
import { truthyKeys } from '@/lib/helpers';

import { VECTOR_COLOR_MAPS } from '@/config/color-maps';
import { TERRESTRIAL_LANDUSE_COLORS } from '@/config/solutions/colors';
import { getSolutionsDataAccessor } from '@/config/solutions/data-access';
import { getTerrestrialDataFormats } from '@/config/solutions/data-formats';
import { LandUseOption, TerrestrialLocationFilterType } from '@/config/solutions/domains';
import { sectionStyleValueState, sectionVisibilityState } from '@/state/sections';
import { landuseFilterState } from '@/state/solutions/landuse-tree';
import {
  TerrestrialFilters,
  TerrestrialLocationFilters,
  terrestrialFiltersState,
} from '@/state/solutions/terrestrial-filters';

export function landuseColorMap(x: string) {
  return TERRESTRIAL_LANDUSE_COLORS[x].css;
}

export const terrestrialColorSpecState = selector<ColorSpec>({
  key: 'terrestrialColorSpecState',
  get: ({ get }) => {
    const style = get(sectionStyleValueState('terrestrial'));

    if (style === 'elevation') {
      return VECTOR_COLOR_MAPS.terrestrialElevation;
    } else if (style === 'slope') {
      return VECTOR_COLOR_MAPS.terrestrialSlope;
    } else {
      // land use will not have a colorSpec, because it's categorical
      return null;
    }
  },
});

export const terrestrialColorFnState = selector<Accessor<string>>({
  key: 'terrestrialColorFnState',
  get: ({ get }) => {
    const style = get(sectionStyleValueState('terrestrial'));

    if (style === 'landuse') {
      return landuseColorMap;
    } else {
      const colorSpec = get(terrestrialColorSpecState);

      if (colorSpec) {
        return colorMap(colorSpec);
      }
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

function landuseFilterValue(p, landuseFilters: Set<LandUseOption>) {
  return landuseFilters.has(p.landuse_desc) ? 1 : 0;
}

function locationFilterValue(p, locationFiltersKeys: TerrestrialLocationFilterType[]) {
  //eslint-disable-next-line eqeqeq -- values are currently sometimes 1 and sometimes true
  return locationFiltersKeys.every((key) => p[key] == true) ? 1 : 0;
}

const landuseFilterSetState = selector<Set<LandUseOption>>({
  key: 'landuseFilterSetState',
  get: ({ get }) => new Set(truthyKeys(get(landuseFilterState))),
});

const locationFilterState = selector<TerrestrialLocationFilters>({
  key: 'locationFilterState',
  get: ({ get }) => get(terrestrialFiltersState).location_filters,
});

const locationFilterKeysState = selector<TerrestrialLocationFilterType[]>({
  key: 'locationFilterKeysState',
  get: ({ get }) => truthyKeys(get(locationFilterState)),
});

function terrestrialFilters(
  filters: TerrestrialFilters,
  landuseFilterSet: Set<LandUseOption>,
  locationFilterKeys: TerrestrialLocationFilterType[],
  doFilter: boolean = true,
) {
  return {
    getFilterValue: ({ properties }) => [
      doFilter ? landuseFilterValue(properties, landuseFilterSet) : 1,
      doFilter ? properties.slope_degrees : 1,
      doFilter ? properties.elevation_m : 1,
      doFilter ? locationFilterValue(properties, locationFilterKeys) : 1,
    ],
    filterRange: [
      [1, 1],
      doFilter ? [...filters.slope_degrees] : [1, 1],
      doFilter ? [...filters.elevation_m] : [1, 1],
      [1, 1],
    ],

    updateTriggers: {
      getFilterValue: [doFilter, landuseFilterSet, locationFilterKeys],
    },

    extensions: [new DataFilterExtension({ filterSize: 4 })],
  };
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

    const colorSpec = get(terrestrialColorSpecState);

    const landuseFilterSet = get(landuseFilterSetState);
    const locationFilterKeys = get(locationFilterKeysState);

    return {
      id: 'terrestrial',
      group: null,
      interactionGroup: 'solutions',
      styleParams: colorSpec && {
        colorMap: {
          fieldSpec,
          colorSpec,
        },
      },
      dataAccessFn: getSolutionsDataAccessor,
      dataFormatsFn: getTerrestrialDataFormats,
      fn: ({ deckProps, zoom, selection }) => {
        const switchoverZoom = 14.5;

        return [
          selectableMvtLayer(
            {
              selectionOptions: {
                selectedFeatureId: selection?.target?.feature.id,
                polygonOffset: -1000,
              },
            },
            deckProps,
            {
              id: `${deckProps.id}@points`,
              data: '/vector/data/natural_terrestrial_combined_points.json',
              visible: zoom < switchoverZoom,
              binary: false,
              filled: true,
              pointAntialiasing: false,
              stroked: false,
            },
            {
              getPointRadius: 35,
              pointRadiusUnit: 'meters',
              pointRadiusMinPixels: 4,
              pointRadiusMaxPixels: 7,
            },
            fillColor(dataColorMap(dataFn, colorFn)),
            terrestrialFilters(filters, landuseFilterSet, locationFilterKeys, zoom < switchoverZoom),
          ),
          selectableMvtLayer(
            {
              selectionOptions: {
                selectedFeatureId: selection?.target?.feature.id,
                polygonOffset: -1000,
              },
            },
            deckProps,
            {
              data: '/vector/data/natural_terrestrial_combined.json',
              minZoom: 14,
              visible: zoom >= switchoverZoom,
              binary: false,
              filled: true,

              getLineWidth: 1,
              lineWidthUnit: 'pixels',
              lineWidthMinPixels: 1,
            },
            border(),
            fillColor(dataColorMap(dataFn, colorFn)),
            terrestrialFilters(filters, landuseFilterSet, locationFilterKeys, zoom >= switchoverZoom),
          ),
        ];
      },
    };
  },
});
