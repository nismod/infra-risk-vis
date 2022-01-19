import { useCallback, useMemo } from 'react';
import _ from 'lodash';

import { LayerDefinition, LayerName, LAYERS } from '../config/layers';
import { ViewName, VIEWS } from '../config/views';
import { DECK_LAYERS } from '../config/deck-layers';
import { boundariesLayer, BoundaryLevel } from '../config/deck-layers/boundaries-layer';
import { selectionLayer } from '../config/deck-layers/selection-layer';
import { labelsLayer } from '../config/deck-layers/labels-layer';

import { VectorHover } from './DataMap';

/**
 * get map style and layers definition based on:
 * - selected background
 * - selected layers, filters
 * - selected data visualisation
 * - any highlights / selections
 */

function getDeckLayersSpec(dataLayerSelection: Record<LayerName, boolean>, view: ViewName) {
  const deckLayers = {};

  for (const layerName of VIEWS[view].layers) {
    if (dataLayerSelection[layerName] == null) continue;

    const layerDefinition = LAYERS[layerName] as LayerDefinition;
    if (layerDefinition == null) throw new Error(`Logical layer '${layerName}' is not defined`);

    const deckLayerSpec = layerDefinition.deckLayer;

    let deckLayerName: string;
    let dataParams: any;
    if (typeof deckLayerSpec === 'object') {
      deckLayerName = deckLayerSpec.baseName;
      dataParams = deckLayerSpec.params;
      if (layerDefinition.getId) {
        deckLayerName = layerDefinition.getId(dataParams);
      }
    } else {
      deckLayerName = deckLayerSpec;
    }

    deckLayers[deckLayerName] = _.merge({}, deckLayers[deckLayerName], {
      visibility: {
        [layerName]: !!dataLayerSelection[layerName],
      },
      params: dataParams ?? {},
      sourceLogicalLayers: [layerName],
    });
  }

  return deckLayers;
}

function getDeckLayers(
  deckLayersSpec: Record<string, any>,
  zoom: number,
  styleParams: any,
  selectedFeature: VectorHover,
  showLabels: boolean,
  showBoundaries: boolean,
  boundaryLevel: BoundaryLevel,
  isRetina: boolean,
) {
  const resLayers = [];

  if (showBoundaries) {
    resLayers.push(boundariesLayer(boundaryLevel));
  }

  for (const [deckLayerName, allParams] of Object.entries(deckLayersSpec)) {
    const deckLayerConfig = DECK_LAYERS[deckLayerName];
    const anyVisible = Object.values(allParams.visibility).some((x) => x);

    let props: any = {
      id: deckLayerName,
      pickable: true,
      visible: anyVisible,
      minZoom: 3,
      maxZoom: 20,
    };

    if (anyVisible) {
      resLayers.push(deckLayerConfig.fn({ props, ...allParams, zoom, styleParams }));
    }
  }

  if (selectedFeature) {
    const { feature } = selectedFeature;

    resLayers.push(selectionLayer(feature, zoom));
  }

  if (showLabels) {
    resLayers.push(labelsLayer(isRetina));
  }

  return resLayers;
}

export function useDeckLayersSpec(dataLayerSelection, view) {
  return useMemo(() => getDeckLayersSpec(dataLayerSelection, view), [dataLayerSelection, view]);
}

export function useMapLayersFunction(
  deckLayersSpec,
  styleParams,
  selectedFeature,
  showLabels,
  showBoundaries,
  boundaryLevel,
  isRetina,
) {
  return useCallback(
    ({ zoom }) =>
      getDeckLayers(
        deckLayersSpec,
        zoom,
        styleParams,
        selectedFeature,
        showLabels,
        showBoundaries,
        boundaryLevel,
        isRetina,
      ),
    [deckLayersSpec, styleParams, selectedFeature, showLabels, showBoundaries, boundaryLevel, isRetina],
  );
}
