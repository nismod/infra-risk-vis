import { MapboxGeoJSONFeature } from 'mapbox-gl';
import { useMemo } from 'react';
import { TileLayer, BitmapLayer, MVTLayer } from 'deck.gl';
import GL from '@luma.gl/constants';

import { LayerDefinition, LayerName, layers } from '../config/layers';
import { BackgroundName } from '../config/types';
import { ViewName, views } from '../config/views';
import { colorToRGB } from '../helpers';

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

// function makeSources<T>(values: T[], keyTransform: (v: T) => string, valueTransform: (v: T) => object) {
//   return Object.fromEntries(values.map((v) => [keyTransform(v), valueTransform(v)]));
// }

function getMapboxSources() {
  const res = {
    satellite: {
      type: 'raster',
      url: 'mapbox://mapbox.satellite',
      tileSize: 256,
    },
    light: {
      type: 'raster',
      tiles: [
        'https://tiles.basemaps.cartocdn.com/light_nolabels/{z}/{x}/{y}.png',
        'https://a.basemaps.cartocdn.com/light_nolabels/{z}/{x}/{y}.png',
        'https://b.basemaps.cartocdn.com/light_nolabels/{z}/{x}/{y}.png',
        'https://c.basemaps.cartocdn.com/light_nolabels/{z}/{x}/{y}.png',
        'https://d.basemaps.cartocdn.com/light_nolabels/{z}/{x}/{y}.png',
      ],
      tileSize: 256,
    },

    // ...makeSources(
    //   ['road_edges', 'bridges', 'elec_edges', 'elec_nodes', 'rail_edges', 'rail_nodes', 'pot_edges', 'abs_nodes'],
    //   (v) => v,
    //   (v) => ({
    //     type: 'vector',
    //     url: `http://localhost:8080/data/${v}.json`,
    //   }),
    // ),
  };

  return res;
}

function visible(isVisible: boolean): 'visible' | 'none' {
  return isVisible ? 'visible' : 'none';
}

function getMapboxLayers({ background, view, dataLayerSelection, highlightedFeature }: MapParams) {
  let res = [];

  res.push({
    id: 'bg-satellite',
    source: 'satellite',
    type: 'raster',
    'source-layer': 'mapbox_satellite_full',
    layout: {
      visibility: visible(background === 'satellite'),
    },
  });

  res.push({
    id: 'bg-light',
    source: 'light',
    type: 'raster',
    layout: {
      visibility: visible(background === 'light'),
    },
  });

  // let highlightLayer: object;
  // let highlightedLayerName: string;
  // if (highlightedFeature) {
  //   const layer = highlightedFeature.layer;

  //   let visualStyle: object;

  //   if (layer.type === 'line') {
  //     visualStyle = {
  //       layout: {
  //         'line-join': 'round',
  //         'line-cap': 'round',
  //       },
  //       paint: {
  //         'line-color': 'yellow',
  //         'line-width': {
  //           base: 1,
  //           stops: [
  //             [3, 1],
  //             [10, 8],
  //             [17, 16],
  //           ],
  //         },
  //       },
  //     };
  //   } else if (layer.type === 'circle') {
  //     visualStyle = {
  //       paint: {
  //         'circle-color': 'yellow',
  //         'circle-radius': {
  //           base: 1,
  //           stops: [
  //             [3, 4],
  //             [10, 12],
  //             [17, 20],
  //           ],
  //         },
  //       },
  //     };
  //   }
  //   if (visualStyle) {
  //     highlightedLayerName = layer.id;
  //     highlightLayer = {
  //       id: `feature-highlight-${highlightedFeature.id}`,
  //       source: layer.source,
  //       'source-layer': layer['source-layer'],
  //       type: layer.type,
  //       filter: ['==', ['id'], highlightedFeature.id],
  //       ...visualStyle,
  //     };
  // }
  // }

  // for (const layerName of views[view].layers) {
  //   const layerDefinition = layers[layerName];

  //   // raster layers are handled by Deck.GL
  //   if (layerDefinition.type === 'raster') continue;

  //   if (layerName === highlightedLayerName) {
  //     res.push(highlightLayer);
  //   }

  //   // TODO: prevent overwriting layout from original style - perform deep merge instead
  //   res.push({
  //     ...layerDefinition.style,
  //     layout: { visibility: visible(dataLayerSelection[layerName]) },
  //   });
  // }

  return res;
}

export function useMapboxStyle(params: MapParams) {
  const sources = useMemo(() => getMapboxSources(), []);
  const layers = useMemo(() => getMapboxLayers(params), [params]);

  const res = useMemo(
    () => ({
      version: 8,
      sources,
      layers,
      glyphs: 'http://localhost:8080/{fontstack}/{range}.pbf',
    }),
    [sources, layers],
  );

  return res;
}

function getDeckLayers(dataLayerSelection: Record<LayerName, boolean>, view: ViewName, onLayerHover) {
  const resLayers = [];

  for (const layerName of views[view].layers) {
    const layerDefinition = layers[layerName] as LayerDefinition;

    if (layerDefinition.type === 'raster') {
      resLayers.push(
        new TileLayer({
          id: layerName,
          data: layerDefinition.sourceUrl,
          visible: !!dataLayerSelection[layerName],
          refinementStrategy: 'no-overlap',
          pickable: true,
          onHover: onLayerHover,
          renderSubLayers: (props) => {
            const {
              bbox: { west, south, east, north },
            } = props.tile;

            return new BitmapLayer(props, {
              data: null,
              image: props.data,
              bounds: [west, south, east, north],
              textureParameters: {
                [GL.TEXTURE_MIN_FILTER]: GL.NEAREST,
                [GL.TEXTURE_MAG_FILTER]: GL.NEAREST,
              },
            });
          },
        }),
      );
    } else if (layerDefinition.type === 'line' || layerDefinition.type === 'circle') {
      const color = colorToRGB(layerDefinition.color);
      let styleProperties: object;

      if (layerDefinition.type === 'line') {
        const { value, unit, minPixels, maxPixels } = layerDefinition.width;
        styleProperties = {
          getLineWidth: value,
          lineWidthUnits: unit,
          lineWidthMinPixels: minPixels,
          lineWidthMaxPixels: maxPixels,
          lineJointRounded: true,
          lineCapRounded: true,
        };
      } else {
        const { value, unit, minPixels, maxPixels } = layerDefinition.radius;
        styleProperties = {
          getPointRadius: value,
          pointRadiusUnits: unit,
          pointRadiusMinPixels: minPixels,
          pointRadiusMaxPixels: maxPixels,
        };
      }
      resLayers.push(
        new MVTLayer({
          id: layerName,
          data: layerDefinition.sourceUrl,
          visible: !!dataLayerSelection[layerName],
          binary: true,
          pickable: true,
          pickingRadius: 20,
          autoHighlight: true,
          highlightColor: [0, 255, 255],
          minZoom: 3,
          maxZoom: 20,
          getLineColor: color,
          getFillColor: color,
          ...styleProperties,
        } as any),
      );
    }
  }

  return resLayers;
}

export function useDeckLayers(params: MapParams, onLayerHover) {
  const { dataLayerSelection, view } = params;

  return useMemo(() => getDeckLayers(dataLayerSelection, view, onLayerHover), [dataLayerSelection, view, onLayerHover]);
}
