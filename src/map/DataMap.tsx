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
import { MapSearch } from './search/MapSearch';
import { LegendContent } from './legend/LegendContent';
import { MapLayerSelection } from './layers/MapLayerSelection';
import { Box } from '@material-ui/core';
import { placeSearchSelectedResultState } from './search/search-state';
import { useRecoilState, useRecoilValue } from 'recoil';
import { backgroundState, showLabelsState, showBoundariesState, boundaryLevelState } from './layers/layers-state';
import { regionHoverState } from './hover/hover-state';

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

export interface RegionHover {
  type: 'region';
  deckLayer: string;
  feature: any;
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

function processRegionHover(layerId, info): RegionHover {
  const { object } = info;

  return {
    type: 'region',
    deckLayer: layerId,
    feature: object,
  };
}

const pickingRadius = 8;
const rasterRegex = /^(coastal|fluvial|surface|cyclone)/;
const boundaryLayerIds = ['boundaries-parish', 'boundaries-community', 'boundaries-enumeration'];

export const DataMap = ({ view, layerSelection, styleParams }) => {
  const background = useRecoilValue(backgroundState);
  const showLabels = useRecoilValue(showLabelsState);
  const showBoundaries = useRecoilValue(showBoundariesState);
  const boundaryLevel = useRecoilValue(boundaryLevelState);

  const [hoveredVectors, setHoveredVectors] = useState<VectorHover[]>([]);
  const [hoveredRasters, setHoveredRasters] = useState<RasterHover[]>([]);
  const [hoveredRegion, setHoveredRegion] = useRecoilState(regionHoverState);
  const [hoverXY, setHoverXY] = useState<[number, number]>(null);

  const [selectedFeature, setSelectedFeature] = useState<VectorHover>(null);

  const deckLayersSpec = useDeckLayersSpec(layerSelection, view);

  const deckIds = useMemo(() => Object.keys(deckLayersSpec), [deckLayersSpec]);
  const rasterLayerIds = useMemo(() => deckIds.filter((l: string) => l.match(rasterRegex)), [deckIds]);
  const vectorLayerIds = useMemo(() => deckIds.filter((l: string) => !l.match(rasterRegex)), [deckIds]);

  // taken from Leaflet source: https://github.com/Leaflet/Leaflet/blob/ee71642691c2c71605bacff69456760cfbc80a2a/src/core/Browser.js#L119
  var isRetina =
    (window.devicePixelRatio || (window.screen as any).deviceXDPI / (window.screen as any).logicalXDPI) > 1;

  const onHover = useCallback(
    (info: any, deck: DeckGL) => {
      const { x, y } = info;

      const newHoveredVectors: VectorHover[] = [];
      const newHoveredRasters: RasterHover[] = [];
      let newPickedRegion: RegionHover = null;

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

      const pickedRegion = deck.pickObject({
        x,
        y,
        layerIds: boundaryLayerIds,
      });
      if (pickedRegion) {
        const layerId = pickedRegion.layer.id;
        newPickedRegion = processRegionHover(layerId, pickedRegion);
      }

      setHoveredVectors(newHoveredVectors);
      setHoveredRasters(newHoveredRasters);
      setHoveredRegion(newPickedRegion);
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

  const deckLayersFunction = useMapLayersFunction(
    deckLayersSpec,
    styleParams,
    selectedFeature,
    showLabels,
    showBoundaries,
    boundaryLevel,
    isRetina,
    background,
  );

  const selectedSearchResult = useRecoilValue(placeSearchSelectedResultState);
  const searchBounds = selectedSearchResult?.boundingBox;

  return (
    <>
      <MapViewport
        layersFunction={deckLayersFunction}
        background={background}
        onHover={onHover}
        onClick={onClick}
        pickingRadius={pickingRadius}
        targetBounds={searchBounds}
      >
        <MapTooltip tooltipXY={hoverXY}>
          {hoveredRasters.length || hoveredVectors.length ? (
            <TooltipContent
              hoveredVectors={hoveredVectors}
              hoveredRasters={hoveredRasters}
              hoveredRegion={hoveredRegion}
            />
          ) : null}
        </MapTooltip>
      </MapViewport>
      <Box position="absolute" top={0} left={0} ml={3} m={1} zIndex={1000}>
        <Box mt={1}>
          <MapSearch />
        </Box>
        <Box mt={1}>
          <MapLayerSelection />
        </Box>
      </Box>
      <Box position="absolute" bottom={0} left={0} m={1} ml={3} zIndex={1000}>
        <MapLegend>
          <LegendContent deckLayersSpec={deckLayersSpec} styleParams={styleParams} />
        </MapLegend>
      </Box>
      {selectedFeature && <FeatureSidebar featureSelection={selectedFeature} />}
    </>
  );
};
