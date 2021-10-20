import { useCallback, useMemo } from 'react';
import _ from 'lodash';

import { LayerDefinition, LayerName, LAYERS } from '../config/layers';
import { ViewName, VIEWS } from '../config/views';
import { DECK_LAYERS } from '../config/deck-layers';

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

const damageMapProps = {
  styleParams: { colorMap: { colorScheme: 'damages', colorField: 'cyclone__rcp_4.5__epoch_2050__conf_None' } },
};

function getDeckLayers(deckLayersSpec: Record<string, any>, zoom: number, showDamages: boolean) {
  const resLayers = [];
  const damageProps = showDamages ? damageMapProps : {};

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

    resLayers.push(deckLayerConfig.fn({ props, ...allParams, zoom, ...damageProps }));
  }

  return resLayers;
}

export function useDeckLayersSpec(dataLayerSelection, view) {
  return useMemo(() => getDeckLayersSpec(dataLayerSelection, view), [dataLayerSelection, view]);
}

export function useMapLayersFunction(deckLayersSpec, showDamages) {
  return useCallback(({ zoom }) => getDeckLayers(deckLayersSpec, zoom, showDamages), [deckLayersSpec, showDamages]);
}
