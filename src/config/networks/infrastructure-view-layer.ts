import { vectorDeckLayer } from 'lib/deck-layers/vector-deck-layer';
import { ViewLayer } from 'lib/data-map/view-layers';
import { ASSETS_SOURCE } from './sources';

export function infrastructureViewLayer(assetId: string, customFn: ({ zoom, styleParams }) => object[]): ViewLayer {
  return {
    id: assetId,
    group: 'infrastructure',
    spatialType: 'vector',
    interactionGroup: 'assets',
    params: {
      assetId,
    },
    fn: ({ deckProps, zoom, styleParams, selection }) =>
      vectorDeckLayer(
        { selectedFeatureId: selection?.target.feature.id, polygonOffset: -1000 },
        deckProps,
        {
          data: ASSETS_SOURCE.getDataUrl({ assetId }),
        },
        ...customFn({ zoom, styleParams }),
      ),
  };
}
