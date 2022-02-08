import { useState } from 'react';
import { useRecoilValue } from 'recoil';
import { Box } from '@mui/material';
import { AttributionControl, NavigationControl, ScaleControl } from 'react-map-gl';

import { MapTooltip } from 'lib/map/MapTooltip';
import { MapViewport } from 'lib/map/MapViewport';
import { MapSearch } from 'lib/map/place-search/MapSearch';
import { placeSearchSelectedResultState } from 'lib/map/place-search/search-state';
import { InteractionGroupConfig, useInteractions } from 'lib/map/interactions/use-interactions';

import { FeatureSidebar } from '../features/FeatureSidebar';
import { useViewLayersSpec, useMapLayersFunction } from './use-map-layers';
import { TooltipContent } from './tooltip/TooltipContent';
import { MapLegend } from './legend/MapLegend';
import { MapLayerSelection } from './layers/MapLayerSelection';
import { LegendContent } from './legend/LegendContent';
import { backgroundState, showLabelsState, showBoundariesState, boundaryLevelState } from './layers/layers-state';
import { useBackgroundConfig } from './use-background-config';
import { INTERACTION_GROUPS } from 'config/interaction-groups';

export const DataMap = ({ view, layerSelection, styleParams }) => {
  const background = useRecoilValue(backgroundState);
  const showLabels = useRecoilValue(showLabelsState);
  const showBoundaries = useRecoilValue(showBoundariesState);
  const boundaryLevel = useRecoilValue(boundaryLevelState);

  const [selectedFeature, setSelectedFeature] = useState(null);

  const viewLayersSpec = useViewLayersSpec(layerSelection, view);

  // taken from Leaflet source: https://github.com/Leaflet/Leaflet/blob/ee71642691c2c71605bacff69456760cfbc80a2a/src/core/Browser.js#L119
  var isRetina =
    (window.devicePixelRatio || (window.screen as any).deviceXDPI / (window.screen as any).logicalXDPI) > 1;

  const { onHover, layerFilter, pickingRadius } = useInteractions(
    viewLayersSpec,
    INTERACTION_GROUPS as Record<string, InteractionGroupConfig>,
  );

  const deckLayersFunction = useMapLayersFunction(
    viewLayersSpec,
    styleParams,
    selectedFeature,
    showLabels,
    showBoundaries,
    boundaryLevel,
    isRetina,
    background,
  );

  const backgroundStyle = useBackgroundConfig(background);

  const selectedSearchResult = useRecoilValue(placeSearchSelectedResultState);
  const searchBounds = selectedSearchResult?.boundingBox;

  return (
    <>
      <MapViewport
        initialViewState={{
          latitude: 18.14,
          longitude: -77.28,
          zoom: 8,
          minZoom: 3,
          maxZoom: 16,
          maxPitch: 0,
        }}
        layersFunction={deckLayersFunction}
        backgroundStyle={backgroundStyle}
        onHover={onHover}
        // onClick={onClick}
        layerRenderFilter={layerFilter}
        pickingRadius={pickingRadius}
        targetBounds={searchBounds}
      >
        <AttributionControl
          customAttribution='Background map data &copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors, style &copy; <a href="https://carto.com/attributions">CARTO</a>. Satellite imagery: <a href="https://s2maps.eu">Sentinel-2 cloudless - https://s2maps.eu</a> by <a href="https://eox.at">EOX IT Services GmbH</a> (Contains modified Copernicus Sentinel data 2020)'
          compact={false}
          style={{
            right: 0,
            bottom: 0,
          }}
        />
        <NavigationControl
          showCompass={false}
          capturePointerMove={true}
          style={{
            right: 10,
            top: 10,
          }}
        />
        <ScaleControl
          maxWidth={100}
          unit="metric"
          style={{
            right: 10,
            bottom: 25,
          }}
        />
        <MapTooltip>
          <TooltipContent />
        </MapTooltip>
      </MapViewport>
      <Box position="absolute" top={0} left={0} ml={3} m={1} zIndex={1000}>
        <Box mb={1}>
          <MapSearch />
        </Box>
        <Box mb={1}>
          <MapLayerSelection />
        </Box>
      </Box>
      <Box position="absolute" bottom={0} left={0} m={1} ml={3} zIndex={1000}>
        <MapLegend>
          <LegendContent viewLayersSpec={viewLayersSpec} styleParams={styleParams} />
        </MapLegend>
      </Box>
      {selectedFeature && <FeatureSidebar featureSelection={selectedFeature} />}
    </>
  );
};
