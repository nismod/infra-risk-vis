import { useRecoilValue } from 'recoil';
import { AttributionControl, NavigationControl, ScaleControl } from 'react-map-gl';

import { DataMap } from 'lib/data-map/DataMap';
import { placeSearchSelectedResultState } from 'lib/map/place-search/search-state';

import { INTERACTION_GROUPS } from 'config/interaction-groups';

import { backgroundState } from './layers/layers-state';
import { useBackgroundConfig } from './use-background-config';
import { MapBoundsFitter } from 'lib/map/MapBoundsFitter';
import { useSaveViewLayers, viewLayersParamsState } from 'state/layers/view-layers-params';
import { viewLayersFlatState } from 'state/layers/view-layers-flat';
import { useEffect } from 'react';
import { DataMapTooltip } from 'lib/data-map/DataMapTooltip';
import { TooltipContent } from './tooltip/TooltipContent';
import { Box } from '@mui/material';
import { MapSearch } from 'lib/map/place-search/MapSearch';
import { MapLayerSelection } from './layers/MapLayerSelection';
import { MapLegend } from './legend/MapLegend';
import { LegendContent } from './legend/LegendContent';
import { FeatureSidebar } from 'features/FeatureSidebar';
import { globalStyleVariables } from 'theme';

export const MapView = ({ view }) => {
  const background = useRecoilValue(backgroundState);
  const viewLayers = useRecoilValue(viewLayersFlatState);
  const saveViewLayers = useSaveViewLayers();

  useEffect(() => {
    saveViewLayers(viewLayers);
  }, [saveViewLayers, viewLayers]);

  const viewLayersParams = useRecoilValue(viewLayersParamsState);

  const backgroundStyle = useBackgroundConfig(background);

  const selectedSearchResult = useRecoilValue(placeSearchSelectedResultState);
  const searchBounds = selectedSearchResult?.boundingBox;

  return (
    <DataMap
      initialViewState={{
        latitude: 18.14,
        longitude: -77.28,
        zoom: 8,
        minZoom: 3,
        maxZoom: 16,
        maxPitch: 0,
      }}
      viewLayers={viewLayers}
      viewLayersParams={viewLayersParams}
      interactionGroups={INTERACTION_GROUPS}
      backgroundStyle={backgroundStyle}
      uiOverlays={
        <>
          <DataMapTooltip>
            <TooltipContent />
          </DataMapTooltip>
          <Box position="absolute" top={0} left={globalStyleVariables.sidebarWidth + 10} ml={3} m={1} zIndex={1000}>
            <Box mb={1}>
              <MapSearch />
            </Box>
            <Box mb={1}>
              <MapLayerSelection />
            </Box>
          </Box>
          <Box position="absolute" bottom={0} left={globalStyleVariables.sidebarWidth + 10} m={1} ml={3} zIndex={1000}>
            <MapLegend>
              <LegendContent />
            </MapLegend>
          </Box>
          <FeatureSidebar />
        </>
      }
    >
      <MapBoundsFitter boundingBox={searchBounds} />

      <AttributionControl
        customAttribution='Background map data &copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors, style &copy; <a href="https://carto.com/attributions">CARTO</a>. Satellite imagery: <a href="https://s2maps.eu">Sentinel-2 cloudless - https://s2maps.eu</a> by <a href="https://eox.at">EOX IT Services GmbH</a> (Contains modified Copernicus Sentinel data 2020)'
        compact={true}
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
          bottom: 40,
        }}
      />
    </DataMap>
  );
};
