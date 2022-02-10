import { infrastructureLayer } from 'config/deck-layers/infrastructure-layer';
import { ViewLayer } from 'lib/view-layers';

export function infrastructureViewLayer(id: string, customFn: ({ zoom, styleParams }) => object[]): ViewLayer {
  return {
    id,
    spatialType: 'vector',
    interactionGroup: 'assets',
    fn: ({ props, zoom, styleParams, selection }) =>
      infrastructureLayer(
        { selectedFeatureId: selection?.target.feature.id },
        props,
        {
          data: `/vector/data/${id}.json`,
        },
        ...customFn({ zoom, styleParams }),
      ),
  };
}
