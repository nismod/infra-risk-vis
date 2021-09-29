import { MapboxGeoJSONFeature } from 'mapbox-gl';
import { useCallback, useMemo, useState } from 'react';

import { readPixelsToArray } from '@luma.gl/core';

import { MapParams, useMapLayersFunction } from './use-map-layers';
import { MapViewport } from './MapViewport';
import { MapTooltip } from './tooltip/MapTooltip';
import { FeatureSidebar } from '../FeatureSidebar';
import { TooltipContent } from './tooltip/TooltipContent';
import DeckGL from 'deck.gl';
import { DECK_LAYERS } from '../config/deck-layers';
import _ from 'lodash';

export interface RasterHover {
  type: 'raster';
  deckLayer: string;
  color: any;
  info: any;
}

export interface VectorHover {
  type: 'vector';
  deckLayer: string;
  feature: any;
  info: any;
}

export type HoveredObject = VectorHover | RasterHover;

export const DataMap = ({ background, view, layerSelection }) => {
  const [hoveredObjects, setHoveredObjects] = useState<HoveredObject[]>([]);
  const [hoverXY, setHoverXY] = useState<[number, number]>(null);

  const [selectedFeatures, setSelectedFeatures] = useState<MapboxGeoJSONFeature[]>([]);

  const mapContentParams = useMemo<MapParams>(
    () => ({
      background,
      view,
      dataLayerSelection: layerSelection,
      highlightedFeature: selectedFeatures?.[0],
    }),
    [background, view, layerSelection, selectedFeatures],
  );

  const onHover = useCallback((info: any, deck: DeckGL) => {
    const { x, y } = info;

    const newHoveredObjects: HoveredObject[] = [];

    if (info.object || info.bitmap) {
      const pickedObjects = deck.pickMultipleObjects({ x, y, radius: 20 });
      for (const picked of pickedObjects) {
        const layerId = picked.layer.id;
        const deckLayerDefinition = DECK_LAYERS[layerId];
        if (deckLayerDefinition.spatialType === 'raster') {
          const { bitmap, sourceLayer } = picked;
          if (bitmap) {
            const pixelColor = readPixelsToArray(sourceLayer.props.image, {
              sourceX: bitmap.pixel[0],
              sourceY: bitmap.pixel[1],
              sourceWidth: 1,
              sourceHeight: 1,
              sourceType: undefined,
            });
            if (pixelColor[3]) {
              newHoveredObjects.push({
                type: 'raster',
                deckLayer: layerId,
                color: pixelColor,
                info: picked,
              });
            }
          }
        } else {
          const { object } = picked;
          newHoveredObjects.push({
            type: 'vector',
            deckLayer: layerId,
            feature: object,
            info: picked,
          });
        }
      }
    }

    setHoveredObjects(newHoveredObjects);
    setHoverXY([x, y]);
  }, []);

  const onLayerClick = useCallback((info: any) => {
    if (!info.bitmap && info.object) {
      setSelectedFeatures([info.object]);
    } else {
      setSelectedFeatures([]);
    }
  }, []);

  const deckLayersFunction = useMapLayersFunction(mapContentParams, onHover);

  return (
    <>
      <MapViewport layersFunction={deckLayersFunction} background={background} onHover={onHover} onClick={onLayerClick}>
        <MapTooltip tooltipXY={hoverXY}>
          {hoveredObjects.length ? <TooltipContent hoveredObjects={hoveredObjects} /> : null}
        </MapTooltip>
      </MapViewport>
      {selectedFeatures.length !== 0 && <FeatureSidebar feature={selectedFeatures[0]} />}
    </>
  );
};
