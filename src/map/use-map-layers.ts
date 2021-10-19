import { MapboxGeoJSONFeature } from 'mapbox-gl';
import { useCallback, useMemo } from 'react';
import _ from 'lodash';

import { LayerDefinition, LayerName, LAYERS } from '../config/layers';
import { BackgroundName } from '../config/backgrounds';
import { ViewName, VIEWS } from '../config/views';
import { DECK_LAYERS } from '../config/deck-layers';

/**
 * get map style and layers definition based on:
 * - selected background
 * - selected layers, filters
 * - selected data visualisation
 * - any highlights / selections
 */

export interface MapParams {
  background: BackgroundName;
  view: ViewName;
  dataLayerSelection: Record<LayerName, boolean>;
  highlightedFeature: MapboxGeoJSONFeature;
}

function getDeckLayersSpec(dataLayerSelection: Record<LayerName, boolean>, view: ViewName) {
  const deckLayers = {};

  for (const layerName of VIEWS[view].layers) {
    if (dataLayerSelection[layerName] == undefined) continue;
    const layerDefinition = LAYERS[layerName] as LayerDefinition;
    if (layerDefinition == undefined) throw new Error(`Logical layer '${layerName}' is not defined`);
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

    deckLayers[deckLayerName] = _.merge(deckLayers[deckLayerName] ?? {}, {
      visibility: {
        [layerName]: !!dataLayerSelection[layerName],
      },
      params: dataParams ?? {},
    });
  }

  return deckLayers;
}

function getDeckLayers(deckLayersSpec: Record<string, any>, zoom: number) {
  const resLayers = [];
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

    resLayers.push(deckLayerConfig.fn({ props, ...allParams, zoom }));
  }

  return resLayers;
}

export function useDeckLayersSpec(dataLayerSelection, view) {
  return useMemo(() => getDeckLayersSpec(dataLayerSelection, view), [dataLayerSelection, view]);
}

export function useMapLayersFunction(deckLayersSpec) {
  return useCallback(({ zoom }) => getDeckLayers(deckLayersSpec, zoom), [deckLayersSpec]);
}
