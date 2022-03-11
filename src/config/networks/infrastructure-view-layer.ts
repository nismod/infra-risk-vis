import { ViewLayer } from 'lib/data-map/view-layers';
import { assetViewLayer } from 'config/assets/asset-view-layer';

export function infrastructureViewLayer(assetId: string, customFn: ({ zoom, styleParams }) => object[]): ViewLayer {
  return assetViewLayer(
    assetId,
    {
      group: 'assets',
      spatialType: 'vector',
      interactionGroup: 'assets',
    },
    -1000,
    customFn,
  );
}
