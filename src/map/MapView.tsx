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

export const MapView = ({ view, children }) => {
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
      uiOverlays={children}
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
