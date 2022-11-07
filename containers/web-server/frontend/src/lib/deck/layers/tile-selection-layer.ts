import { featureFilter } from '../props/feature-filter';
import { GetColor, fillColor, strokeColor } from '../props/style';
import { geoJsonLayer } from './base';

export interface TileSelectionLayerOptions {
  selectedFeatureId: number | null;
  uniqueIdProperty?: string;
  selectionFillColor?: GetColor;
  selectionLineColor?: GetColor;
  polygonOffset?: number;
}
export function tileSelectionLayer(
  tileProps,
  {
    selectedFeatureId,
    uniqueIdProperty,
    selectionFillColor = [0, 255, 255, 255],
    selectionLineColor = [0, 255, 255, 255],
    polygonOffset = 0,
  }: TileSelectionLayerOptions,
) {
  return geoJsonLayer(
    tileProps,
    {
      id: tileProps.id + '-selection',
      pickable: false,
      getPolygonOffset: ({ layerIndex }) => [0, -layerIndex * 100 + polygonOffset],
      visible: selectedFeatureId != null,
      refinementStrategy: 'no-overlap',

      getLineWidth: 2,
      lineWidthUnits: 'pixels',
    },
    // use on-GPU filter extension to only show the selected feature
    featureFilter(selectedFeatureId, uniqueIdProperty),

    fillColor(selectionFillColor),
    strokeColor(selectionLineColor),
  );
}
