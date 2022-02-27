import { ViewLayer } from 'lib/data-map/view-layers';
import { infrastructureDeckLayer } from 'config/networks/infrastructure-deck-layer';
import { border, vectorColor } from 'lib/deck-layers/utils';
import { RegionLevel } from './metadata';

export function populationViewLayer(boundaryLevel: RegionLevel): ViewLayer {
  return {
    id: `population_${boundaryLevel}`,
    interactionGroup: 'regions',
    spatialType: 'vector',
    group: 'regions',
    params: {
      boundaryLevel,
    },
    fn: ({ deckProps, zoom, styleParams, selection }) =>
      infrastructureDeckLayer(
        { selectedFeatureId: selection?.target.feature.id },
        deckProps,
        {
          data: `/vector/data/regions_${boundaryLevel}.json`,
        },
        zoom > 12 ? border([40, 40, 40, 255]) : {},
        vectorColor('fill', '#ccc', {
          colorMap: { colorScheme: 'population', colorField: 'population_density_per_km2' },
        }),
        {
          highlightColor: [0, 255, 255, 100],
        },
      ),
  };
}
