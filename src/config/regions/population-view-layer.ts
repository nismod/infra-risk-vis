import { ViewLayer } from 'lib/data-map/view-layers';
import { MVTLayer } from 'deck.gl';
import { infrastructureDeckLayer } from 'config/networks/infrastructure-deck-layer';
import { vectorColor } from 'lib/deck-layers/utils';

export function populationViewLayer(): ViewLayer {
  return {
    id: 'population',
    interactionGroup: 'regions',
    spatialType: 'vector',
    group: 'regions',
    fn: ({ deckProps, zoom, styleParams, selection }) =>
      infrastructureDeckLayer(
        { selectedFeatureId: selection?.target.feature.id },
        deckProps,
        {
          data: `/vector/data/population.json`,
        },
        vectorColor('fill', '#ccc', { colorMap: { colorScheme: 'population', colorField: 'population_density_per_km2' } }),
        // ...customFn({ zoom, styleParams }),
      ),
  };
}
