import { ViewLayer, ViewLayerFunctionOptions } from 'lib/data-map/view-layers';
import { assetViewLayer } from 'config/assets/asset-view-layer';
import { assetDataAccessFunction } from 'config/assets/data-accessor';

export function infrastructureViewLayer(
  assetId: string,
  customFn: ({ zoom, styleParams }: ViewLayerFunctionOptions) => object[],
): ViewLayer {
  return assetViewLayer(
    assetId,
    {
      group: 'assets',
      spatialType: 'vector',
      interactionGroup: 'assets',
    },
    -1000,
    customFn,
    assetDataAccessFunction(assetId),
  );
}
