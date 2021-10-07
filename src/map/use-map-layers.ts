import { MapboxGeoJSONFeature } from 'mapbox-gl';
import { useMemo } from 'react';
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

function getDeckLayers(dataLayerSelection: Record<LayerName, boolean>, view: ViewName, onLayerHover, zoom) {
  const deckLayerNames = new Set<string>();
  const deckLayerParams = {};

  for (const layerName of VIEWS[view].layers) {
    const layerDefinition = LAYERS[layerName] as LayerDefinition;
    const deckLayerSpec = layerDefinition.deckLayer;
    let deckLayerName: string;
    let deckLayerDataParams: any;
    if (typeof deckLayerSpec === 'object') {
      deckLayerName = deckLayerSpec.baseName;
      deckLayerDataParams = deckLayerSpec.params;
      if (layerDefinition.getId) {
        deckLayerName = layerDefinition.getId(deckLayerDataParams);
      }
    } else {
      deckLayerName = deckLayerSpec;
    }

    deckLayerNames.add(deckLayerName);
    deckLayerParams[deckLayerName] = _.merge(deckLayerParams[deckLayerName] ?? {}, {
      visibility: {
        [layerName]: !!dataLayerSelection[layerName],
      },
      params: deckLayerDataParams ?? {},
    });
  }

  const resLayers = [];
  for (const deckLayerName of deckLayerNames.values()) {
    const deckLayerConfig = DECK_LAYERS[deckLayerName];
    const allParams = deckLayerParams[deckLayerName];
    const anyVisible = Object.values(allParams.visibility).some((x) => x);

    let props: any = {
      id: deckLayerName,
      pickable: true,
      // onHover: onLayerHover,
      visible: anyVisible,
      minZoom: 3,
      maxZoom: 20,
    };

    if (deckLayerConfig.type === 'MVTLayer') {
      props = { ...props, binary: true, autoHighlight: true, highlightColor: [0, 255, 255] };
    } else if (deckLayerConfig.type === 'TileLayer') {
      // do nothing for now
    }

    resLayers.push(deckLayerConfig.fn({ props, ...deckLayerParams[deckLayerName], zoom }));
  }
  // if (layerDefinition.type === 'raster') {
  //   resLayers.push(
  //     new TileLayer({
  //       id: layerName,
  //       visible: !!dataLayerSelection[layerName],
  //       pickable: true,
  //       onHover: onLayerHover,
  //     }),
  //   );
  // } else if (layerDefinition.type === 'line' || layerDefinition.type === 'circle') {

  //   resLayers.push(
  //     new MVTLayer({
  //       id: layerName,
  //       visible: !!dataLayerSelection[layerName],
  //       binary: true,
  //       pickable: true,
  //       autoHighlight: true,
  //       highlightColor: [0, 255, 255],
  //       minZoom: 3,
  //       maxZoom: 20,
  //     } as any),
  //   );
  // }
  // }

  return resLayers;
}

export function useMapLayersFunction(params: MapParams, onLayerHover) {
  const { dataLayerSelection, view } = params;

  return useMemo(() => {
    return ({ zoom }) => getDeckLayers(dataLayerSelection, view, onLayerHover, zoom);
  }, [dataLayerSelection, view, onLayerHover]);
}
