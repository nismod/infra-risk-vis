import { useCallback, useMemo } from 'react';
import _ from 'lodash';

import { LayerDefinition, LayerName, LAYERS } from '../config/layers';
import { ViewName, VIEWS } from '../config/views';
import { getDeckLayers } from 'config/get-deck-layers';

/**
 * get map style and layers definition based on:
 * - selected background
 * - selected layers, filters
 * - selected data visualisation
 * - any highlights / selections
 */

function getSingleViewLayerConfig({ getId, viewLayer }: LayerDefinition) {
  let name: string;
  let dataParams: any;
  if (typeof viewLayer === 'object') {
    dataParams = viewLayer.params;
    name = getId?.(dataParams) ?? viewLayer.baseName;
  } else {
    name = viewLayer;
  }

  return {
    name,
    dataParams,
  };
}

function getViewLayersSpec(logicalLayerSelection: Record<LayerName, boolean>, view: ViewName) {
  const viewLayerSpec = {};

  for (const logicalLayerName of VIEWS[view].layers) {
    if (logicalLayerSelection[logicalLayerName] == null) continue;

    const logicalLayerDefinition = LAYERS[logicalLayerName] as LayerDefinition;
    if (logicalLayerDefinition == null) throw new Error(`Logical layer '${logicalLayerName}' is not defined`);

    const { name, dataParams } = getSingleViewLayerConfig(logicalLayerDefinition);

    const singleViewLayerSpec = {
      visibility: {
        [logicalLayerName]: !!logicalLayerSelection[logicalLayerName],
      },
      params: dataParams ?? {},
      sourceLogicalLayers: [logicalLayerName],
    };

    viewLayerSpec[name] = _.merge({}, viewLayerSpec[name], singleViewLayerSpec);
  }

  return viewLayerSpec;
}

export function useViewLayersSpec(dataLayerSelection, view) {
  return useMemo(() => getViewLayersSpec(dataLayerSelection, view), [dataLayerSelection, view]);
}

export function useMapLayersFunction(
  viewLayersSpec,
  styleParams,
  selectedFeature,
  showLabels,
  showBoundaries,
  boundaryLevel,
  isRetina,
  background,
) {
  return useCallback(
    ({ zoom }) =>
      getDeckLayers(
        viewLayersSpec,
        zoom,
        styleParams,
        selectedFeature,
        showLabels,
        showBoundaries,
        boundaryLevel,
        isRetina,
        background,
      ),
    [viewLayersSpec, styleParams, selectedFeature, showLabels, showBoundaries, boundaryLevel, isRetina, background],
  );
}
