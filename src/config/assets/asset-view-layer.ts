import {
  DataAccess,
  StyleParams,
  ViewLayer,
  ViewLayerDataFunctionOptions,
  ViewLayerFunctionOptions,
} from 'lib/data-map/view-layers';
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
  customFn: ({ zoom, styleParams }: { zoom: number; styleParams?: StyleParams }) => object[],
  customDataAccessFn: ({ styleParams }: ViewLayerDataFunctionOptions) => DataAccess,
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
    fn: ({ deckProps, zoom, styleParams, selection }: ViewLayerFunctionOptions) =>
      selectableMvtLayer(
        {
          selectionOptions: {
            selectedFeatureId: selection?.target.feature.id,
            polygonOffset: selectionPolygonOffset,
          },
          dataLoaderOptions: {
            dataLoader: customDataAccessFn?.({ styleParams })?.dataLoader,
          },
        },
        deckProps,
        {
          data: ASSETS_SOURCE.getDataUrl({ assetId }),
        },
        ...customFn({ zoom, styleParams }),
      ),
    dataAccessFn: customDataAccessFn,
  };
}
