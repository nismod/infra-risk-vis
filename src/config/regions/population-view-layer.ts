import { ViewLayer } from 'lib/data-map/view-layers';
import { MVTLayer } from 'deck.gl';
import { infrastructureDeckLayer } from 'config/networks/infrastructure-deck-layer';
import { border, vectorColor } from 'lib/deck-layers/utils';
import { BoundaryLevel } from './metadata';

export function populationViewLayer(boundaryLevel: BoundaryLevel): ViewLayer {
  return {
    id: `population_${boundaryLevel}`,
    interactionGroup: 'regions',
    spatialType: 'vector',
    group: 'regions',
    params: {
      boundaryLevel: 'enumeration',
    },
    fn: ({ deckProps, zoom, styleParams, selection }) =>
      infrastructureDeckLayer(
        { selectedFeatureId: selection?.target.feature.id },
        deckProps,
        {
          data: `/vector/data/population.json`,
        },
        border([20, 20, 20, 255]),
        vectorColor('fill', '#ccc', {
          colorMap: { colorScheme: 'population', colorField: 'population_density_per_km2' },
        }),
        {
          highlightColor: [0, 255, 255, 100],
        },
      ),
  };
}
