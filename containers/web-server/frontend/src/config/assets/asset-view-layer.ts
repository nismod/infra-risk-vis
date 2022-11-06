import { colorMap } from '@/lib/color-map';
import {
  StyleParams,
  ViewLayer,
  ViewLayerDataAccessFunction,
  ViewLayerRenderDetailsFunction,
  ViewLayerRenderTooltipFunction,
} from '@/lib/data-map/view-layers';
import { selectableMvtLayer } from '@/lib/deck/layers/selectable-mvt-layer';
import { dataColorMap } from '@/lib/deck/props/color-map';
import { GetColor } from '@/lib/deck/props/style';

import { SOURCES } from '../sources';
import { getAssetDataAccessor } from './data-access';
import { getAssetDataFormats } from './data-formats';

interface ViewLayerMetadata {
  spatialType: string;
  interactionGroup: string;
}

export interface DataStyle {
  getColor?: GetColor;
}
export interface AssetViewLayerCustomFunctionOptions {
  zoom: number;
  dataStyle?: DataStyle;
}
export type AssetViewLayerCustomFunction = (options: AssetViewLayerCustomFunctionOptions) => object[];
export interface AssetViewLayerOptions {
  assetId: string;
  metadata: ViewLayerMetadata;
  selectionPolygonOffset?: number;
  styleParams?: StyleParams;
  customFn?: AssetViewLayerCustomFunction;
  customDataAccessFn?: ViewLayerDataAccessFunction;
  renderTooltip?: ViewLayerRenderTooltipFunction;
  renderDetails?: ViewLayerRenderDetailsFunction;
}

export function assetViewLayer({
  assetId,
  metadata: { spatialType, interactionGroup },
  selectionPolygonOffset = -1000,
  styleParams,
  customFn,
  customDataAccessFn,
  renderTooltip,
  renderDetails,
}: AssetViewLayerOptions): ViewLayer {
  const dataStyle: DataStyle = styleParams?.colorMap
    ? {
        getColor: dataColorMap(
          getAssetDataAccessor(assetId, styleParams.colorMap.fieldSpec),
          colorMap(styleParams.colorMap.colorSpec),
        ),
      }
    : null;

  const dataLoader = customDataAccessFn?.(styleParams?.colorMap?.fieldSpec)?.dataLoader;

  return {
    id: assetId,
    spatialType,
    interactionGroup,
    params: {
      assetId,
    },
    styleParams,
    fn: ({ deckProps, zoom, selection }) =>
      selectableMvtLayer(
        {
          selectionOptions: {
            selectedFeatureId: selection?.target.feature.id,
            polygonOffset: selectionPolygonOffset,
          },
          dataLoaderOptions: {
            dataLoader,
          },
        },
        deckProps,
        {
          data: SOURCES.vector.getUrl(assetId),
        },
        ...(customFn?.({ zoom, dataStyle }) ?? []),
      ),
    dataAccessFn: customDataAccessFn,
    dataFormatsFn: getAssetDataFormats,
    renderTooltip,
    renderDetails,
  };
}
