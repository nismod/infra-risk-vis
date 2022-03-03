import { ViewLayer } from 'lib/data-map/view-layers';
import { vectorDeckLayer } from 'lib/deck-layers/vector-deck-layer';
import { border, vectorColor } from 'lib/deck-layers/utils';
import { RegionLevel } from './metadata';
import { REGIONS_SOURCE } from './source';

export function populationViewLayer(regionLevel: RegionLevel): ViewLayer {
  const source = REGIONS_SOURCE;

  return {
    id: `population_${regionLevel}`,
    interactionGroup: 'regions',
    spatialType: 'vector',
    group: 'regions',
    params: {
      regionLevel,
    },
    fn: ({ deckProps, zoom, styleParams, selection }) =>
      vectorDeckLayer(
        { selectedFeatureId: selection?.target.feature.id },
        deckProps,
        {
          data: source.getDataUrl({ regionLevel }),
        },
        regionLevel === 'parish' || zoom > 12 ? border([40, 40, 40, 255]) : {},
        vectorColor('fill', '#ccc', {
          colorMap: { colorScheme: 'population', colorField: 'population_density_per_km2' },
        }),
        {
          highlightColor: [0, 255, 255, 100],
        },
      ),
  };
}
