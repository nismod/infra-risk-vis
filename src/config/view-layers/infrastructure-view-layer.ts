import { infrastructureLayer } from 'config/deck-layers/infrastructure-layer';

export const infrastructureViewLayer = (id: string, customFn: ({ zoom, styleParams }) => object[]) => {
  return {
    id,
    spatialType: 'vector',
    interactionGroup: 'assets',
    fn: ({ props, zoom, styleParams, selectedFeatureId }) =>
      infrastructureLayer(
        { selectedFeatureId },
        props,
        {
          data: `/vector/data/${id}.json`,
        },
        ...customFn({ zoom, styleParams }),
      ),
  };
};
