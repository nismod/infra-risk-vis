import { MapboxGeoJSONFeature } from 'mapbox-gl';
import { useMemo } from 'react';
import { TileLayer, BitmapLayer, MVTLayer } from 'deck.gl';
import GL from '@luma.gl/constants';

import { LayerDefinition, LayerName, layers } from '../config/layers';
import { BackgroundName } from '../config/backgrounds';
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
