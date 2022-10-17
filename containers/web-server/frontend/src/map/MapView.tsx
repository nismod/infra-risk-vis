import { Box } from '@mui/material';
import { useCallback, useEffect } from 'react';
import { AttributionControl, NavigationControl, ScaleControl } from 'react-map-gl';
import { atom, useRecoilState, useRecoilValue } from 'recoil';
import { globalStyleVariables } from 'theme';

import { BoundingBox } from '@/lib/bounding-box';
import { DataMap } from '@/lib/data-map/DataMap';
import { DataMapTooltip } from '@/lib/data-map/DataMapTooltip';
import { MapBoundsFitter } from '@/lib/map/MapBoundsFitter';
import { MapSearch } from '@/lib/map/place-search/MapSearch';
import { PlaceSearchResult } from '@/lib/map/place-search/use-place-search';

import { interactionGroupsState } from '@/state/layers/interaction-groups';
import { viewLayersFlatState } from '@/state/layers/view-layers-flat';
import { useSaveViewLayers, viewLayersParamsState } from '@/state/layers/view-layers-params';

import { MapLayerSelection } from './layers/MapLayerSelection';
import { backgroundState } from './layers/layers-state';
import { MapLegend } from './legend/MapLegend';
import { TooltipContent } from './tooltip/TooltipContent';
import { useBackgroundConfig } from './use-background-config';

export const mapFitBoundsState = atom<BoundingBox>({
  key: 'mapFitBoundsState',
  default: null,
});

export const MapView = () => {
  const background = useRecoilValue(backgroundState);
  const viewLayers = useRecoilValue(viewLayersFlatState);
  const saveViewLayers = useSaveViewLayers();

  useEffect(() => {
    saveViewLayers(viewLayers);
  }, [saveViewLayers, viewLayers]);

  const viewLayersParams = useRecoilValue(viewLayersParamsState);

  const interactionGroups = useRecoilValue(interactionGroupsState);

  const backgroundStyle = useBackgroundConfig(background);

  const [fitBounds, setFitBounds] = useRecoilState(mapFitBoundsState);
  const handleSelectedSearchResult = useCallback(
    (result: PlaceSearchResult) => {
      setFitBounds(result.boundingBox);
    },
    [setFitBounds],
  );

  return (
    <DataMap
      initialViewState={{
        latitude: 20.0,
        longitude: -40.0,
        zoom: 3,
        minZoom: 2,
        maxZoom: 8,
        maxPitch: 0,
      }}
      viewLayers={viewLayers}
      viewLayersParams={viewLayersParams}
      interactionGroups={interactionGroups}
      backgroundStyle={backgroundStyle}
      uiOverlays={
        <>
          <DataMapTooltip>
            <TooltipContent />
          </DataMapTooltip>
          <Box position="absolute" top={0} left={globalStyleVariables.controlSidebarWidth} ml={3} m={1} zIndex={1000}>
            <Box mb={1}>
              <MapSearch onSelectedResult={handleSelectedSearchResult} />
            </Box>
            <Box mb={1}>
              <MapLayerSelection />
            </Box>
          </Box>
          <Box
            position="absolute"
            bottom={0}
            left={globalStyleVariables.controlSidebarWidth}
            m={1}
            ml={1}
            zIndex={1000}
          >
            <MapLegend />
          </Box>
        </>
      }
    >
      <MapBoundsFitter boundingBox={fitBounds} />

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
