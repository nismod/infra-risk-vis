import { ViewLayer } from 'lib/data-map/view-layers';
import { selectableMvtLayer } from 'lib/deck/layers/selectable-mvt-layer';
import { ASSETS_SOURCE } from './source';

interface ViewLayerMetadata {
  group: string;
  spatialType: string;
  interactionGroup: string;
}

export function assetViewLayer(
  assetId: string,
  metadata: ViewLayerMetadata,
  selectionPolygonOffset: number,
  customFn: ({ zoom, styleParams }) => object[],
): ViewLayer {
  const { group, spatialType, interactionGroup } = metadata;

  return {
    id: assetId,
    group,
    spatialType,
    interactionGroup,
    params: {
      assetId,
    },
    fn: ({ deckProps, zoom, styleParams, selection }) =>
      selectableMvtLayer(
        { selectedFeatureId: selection?.target.feature.id, polygonOffset: selectionPolygonOffset },
        deckProps,
        {
          data: ASSETS_SOURCE.getDataUrl({ assetId }),
        },
        ...customFn({ zoom, styleParams }),
      ),
  };
}
