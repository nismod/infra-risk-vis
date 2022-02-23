import { infrastructureDeckLayer } from 'config/networks/infrastructure-deck-layer';
import { ViewLayer } from 'lib/data-map/view-layers';

export function infrastructureViewLayer(id: string, customFn: ({ zoom, styleParams }) => object[]): ViewLayer {
  return {
    id,
    group: 'infrastructure',
    spatialType: 'vector',
    interactionGroup: 'assets',
    fn: ({ deckProps, zoom, styleParams, selection }) =>
      infrastructureDeckLayer(
        { selectedFeatureId: selection?.target.feature.id },
        deckProps,
        {
          data: `/vector/data/${id}.json`,
        },
        ...customFn({ zoom, styleParams }),
      ),
  };
}
