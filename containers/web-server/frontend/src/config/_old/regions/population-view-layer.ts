import { colorMap } from '@/lib/color-map';
import { FieldSpec, ViewLayer } from '@/lib/data-map/view-layers';
import { selectableMvtLayer } from '@/lib/deck/layers/selectable-mvt-layer';
import { dataColorMap } from '@/lib/deck/props/color-map';
import { featureProperty } from '@/lib/deck/props/data-source';
import { border, fillColor } from '@/lib/deck/props/style';

import { VECTOR_COLOR_MAPS } from '@/config/color-maps';

import { RegionLevel } from './metadata';
import { REGIONS_SOURCE } from './source';

export function populationViewLayer(regionLevel: RegionLevel): ViewLayer {
  const source = REGIONS_SOURCE;

  const fieldSpec: FieldSpec = {
    fieldGroup: 'properties',
    field: 'population_density_per_km2',
  };
  const colorSpec = VECTOR_COLOR_MAPS.population;

  return {
    id: `population_${regionLevel}`,
    interactionGroup: 'regions',
    spatialType: 'vector',
    params: {
      regionLevel,
    },
    styleParams: {
      colorMap: {
        fieldSpec,
        colorSpec,
      },
    },
    dataFormatsFn: () => ({
      getDataLabel: () => 'Population density',
      getValueFormatted: (value: number) => `${value.toLocaleString()}/km²`,
    }),
    fn: ({ deckProps, zoom, selection }) =>
      selectableMvtLayer(
        { selectionOptions: { selectedFeatureId: selection?.target.feature.id } },
        deckProps,
        {
          data: source.getDataUrl({ regionLevel }),
        },
        (regionLevel === 'parish' || zoom > 12) && border([40, 40, 40, 255]),
        fillColor(dataColorMap(featureProperty('population_density_per_km2'), colorMap(colorSpec))),
        {
          highlightColor: [0, 255, 255, 100],
        },
      ),
  };
}
