import { colorMap } from '@/lib/color-map';
import { d3Scale, d3ScaleChromatic, invertColorScale } from '@/lib/data-map/color-maps';
import { ColorSpec, FieldSpec, ViewLayer } from '@/lib/data-map/view-layers';
import { selectableMvtLayer } from '@/lib/deck/layers/selectable-mvt-layer';
import { dataColorMap } from '@/lib/deck/props/color-map';
import { featureProperty } from '@/lib/deck/props/data-source';
import { border, fillColor } from '@/lib/deck/props/style';

import { RegionLevel } from './metadata';
import { REGIONS_SOURCE } from './source';

const POPULATION_COLOR_SPEC: ColorSpec = {
  scale: d3Scale.scaleSequential,
  scheme: invertColorScale(d3ScaleChromatic.interpolateInferno),
  range: [0, 1e4],
  empty: '#ccc',
};

export function populationViewLayer(regionLevel: RegionLevel): ViewLayer {
  const source = REGIONS_SOURCE;

  const fieldSpec: FieldSpec = {
    fieldGroup: 'properties',
    field: 'population_density_per_km2',
  };
  const colorSpec = POPULATION_COLOR_SPEC;

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
