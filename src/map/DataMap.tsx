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
import { VIEWS } from '../config/views';

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

function processRasterHover(layerId, info): RasterHover {
  const { bitmap, sourceLayer } = info;
  if (bitmap) {
    const pixelColor = readPixelsToArray(sourceLayer.props.image, {
      sourceX: bitmap.pixel[0],
      sourceY: bitmap.pixel[1],
      sourceWidth: 1,
      sourceHeight: 1,
      sourceType: undefined,
    });
    if (pixelColor[3]) {
      return {
        type: 'raster',
        deckLayer: layerId,
        color: pixelColor,
        info,
      };
    } else return null;
  }
}

function processVectorHover(layerId, info): VectorHover {
  const { object } = info;

  return {
    type: 'vector',
    deckLayer: layerId,
    feature: object,
    info,
  };
}

export const DataMap = ({ background, view, layerSelection }) => {
  // const [hoveredObjects, setHoveredObjects] = useState<HoveredObject[]>([]);
  const [hoveredVectors, setHoveredVectors] = useState<VectorHover[]>([]);
  const [hoveredRasters, setHoveredRasters] = useState<RasterHover[]>([]);
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

  const viewLayers = useMemo(() => VIEWS[view].layers, [view]);

  const rasterLayerIds = useMemo(() => {
    return viewLayers.filter((l: string) => l.match(/^(coastal|fluvial|surface|cyclone)/));
  }, [viewLayers]);

  const vectorLayerIds = useMemo(() => {
    return viewLayers.filter((l: string) => !l.match(/^(coastal|fluvial|surface|cyclone)/));
  }, [viewLayers]);

  const onHover = useCallback(
    (info: any, deck: DeckGL) => {
      const { x, y } = info;

      const newHoveredVectors: VectorHover[] = [];
      const newHoveredRasters: RasterHover[] = [];
      // if (info.object || info.bitmap) {

      // }

      const pickedVectorInfo = deck.pickObject({ x, y, layerIds: vectorLayerIds, radius: 8 });

      if (pickedVectorInfo) {
        const layerId = pickedVectorInfo.layer.id;
        newHoveredVectors.push(processVectorHover(layerId, pickedVectorInfo));
      }

      const pickedObjects = deck.pickMultipleObjects({ x, y, layerIds: rasterLayerIds });
      for (const picked of pickedObjects) {
        const layerId = picked.layer.id;
        const deckLayerDefinition = DECK_LAYERS[layerId];
        if (deckLayerDefinition.spatialType === 'raster') {
          const rasterHover = processRasterHover(layerId, picked);

          if (rasterHover) newHoveredRasters.push(rasterHover);
        } else {
          // it's not expected to get any non-raster layers here
        }
      }
      setHoveredVectors(newHoveredVectors);
      setHoveredRasters(newHoveredRasters);
      // setHoveredObjects(newHoveredObjects);
      setHoverXY([x, y]);
    },
    [rasterLayerIds, vectorLayerIds],
  );

  const onClick = useCallback(
    (info: any, deck: DeckGL) => {
      const { x, y } = info;
      console.log(deck);
      const pickedVectorInfo = deck.pickObject({ x, y, layerIds: vectorLayerIds, radius: 8 });

      if (pickedVectorInfo) {
        // const layerId = pickedVectorInfo.layer.id;
        setSelectedFeatures([pickedVectorInfo.object]);
        // newHoveredVectors.push(processVectorHover(layerId, pickedVectorInfo));
      } else {
        setSelectedFeatures([]);
      }

      // if (!info.bitmap && info.object) {
      //   setSelectedFeatures([info.object]);
      // } else {
      //   setSelectedFeatures([]);
      // }
    },
    [vectorLayerIds],
  );

  const deckLayersFunction = useMapLayersFunction(mapContentParams);

  return (
    <>
      <MapViewport layersFunction={deckLayersFunction} background={background} onHover={onHover} onClick={onClick}>
        <MapTooltip tooltipXY={hoverXY}>
          {hoveredRasters.length || hoveredVectors.length ? (
            <TooltipContent hoveredVectors={hoveredVectors} hoveredRasters={hoveredRasters} />
          ) : null}
        </MapTooltip>
      </MapViewport>
      {selectedFeatures.length !== 0 && <FeatureSidebar feature={selectedFeatures[0]} />}
    </>
  );
};
