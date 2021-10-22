import { MapboxGeoJSONFeature } from 'mapbox-gl';
import { useCallback, useMemo, useState } from 'react';

import { readPixelsToArray } from '@luma.gl/core';

import { useDeckLayersSpec, useMapLayersFunction } from './use-map-layers';
import { MapViewport } from './MapViewport';
import { MapTooltip } from './tooltip/MapTooltip';
import { FeatureSidebar } from '../features/FeatureSidebar';
import { TooltipContent } from './tooltip/TooltipContent';
import DeckGL from 'deck.gl';
import { DECK_LAYERS } from '../config/deck-layers';
import { MapLegend } from './legend/MapLegend';
import { LegendContent } from './legend/LegendContent';
import { MapLayerSelection } from './MapLayerSelection';

export interface RasterHover {
  type: 'raster';
  deckLayer: string;
  logicalLayer: string;
  color: any;
  info: any;
}

export interface VectorHover {
  type: 'vector';
  deckLayer: string;
  logicalLayer: string;
  feature: any;
  info: any;
}

export interface FeatureSelection {
  feature: MapboxGeoJSONFeature;
  sourceDeckLayer: string;
  sourceLogicalLayer: string;
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
        logicalLayer: DECK_LAYERS[layerId].getLogicalLayer?.({ deckLayerId: layerId, feature: null }) ?? layerId,
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
    logicalLayer: DECK_LAYERS[layerId].getLogicalLayer?.({ deckLayerId: layerId, feature: object }) ?? layerId,
    feature: object,
    info,
  };
}

const pickingRadius = 8;
const rasterRegex = /^(coastal|fluvial|surface|cyclone)/;

export const DataMap = ({ background, view, layerSelection, styleParams, onBackground }) => {
  const [hoveredVectors, setHoveredVectors] = useState<VectorHover[]>([]);
  const [hoveredRasters, setHoveredRasters] = useState<RasterHover[]>([]);
  const [hoverXY, setHoverXY] = useState<[number, number]>(null);

  const [selectedFeature, setSelectedFeature] = useState<VectorHover>(null);

  const deckLayersSpec = useDeckLayersSpec(layerSelection, view);

  const deckIds = useMemo(() => Object.keys(deckLayersSpec), [deckLayersSpec]);
  const rasterLayerIds = useMemo(() => deckIds.filter((l: string) => l.match(rasterRegex)), [deckIds]);
  const vectorLayerIds = useMemo(() => deckIds.filter((l: string) => !l.match(rasterRegex)), [deckIds]);

  const onHover = useCallback(
    (info: any, deck: DeckGL) => {
      const { x, y } = info;

      const newHoveredVectors: VectorHover[] = [];
      const newHoveredRasters: RasterHover[] = [];

      const pickedVectorInfo = deck.pickObject({ x, y, layerIds: vectorLayerIds, radius: pickingRadius });

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
      setHoverXY([x, y]);
    },
    [rasterLayerIds, vectorLayerIds],
  );

  const onClick = useCallback(
    (info: any, deck: DeckGL) => {
      const { x, y } = info;
      const pickedVectorInfo = deck.pickObject({ x, y, layerIds: vectorLayerIds, radius: pickingRadius });

      if (pickedVectorInfo) {
        setSelectedFeature(processVectorHover(pickedVectorInfo.layer.id, pickedVectorInfo));
      } else {
        setSelectedFeature(null);
      }
    },
    [vectorLayerIds],
  );

  const deckLayersFunction = useMapLayersFunction(deckLayersSpec, styleParams);

  return (
    <>
      <MapViewport
        layersFunction={deckLayersFunction}
        background={background}
        onHover={onHover}
        onClick={onClick}
        pickingRadius={pickingRadius}
      >
        <MapTooltip tooltipXY={hoverXY}>
          {hoveredRasters.length || hoveredVectors.length ? (
            <TooltipContent hoveredVectors={hoveredVectors} hoveredRasters={hoveredRasters} />
          ) : null}
        </MapTooltip>
      </MapViewport>
      <MapLayerSelection background={background} onBackground={onBackground} />
      {selectedFeature && <FeatureSidebar featureSelection={selectedFeature} />}
      <MapLegend>
        <LegendContent deckLayersSpec={deckLayersSpec} styleParams={styleParams}/>
      </MapLegend>
    </>
  );
};
