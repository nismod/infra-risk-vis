import React, { FC, useCallback, useEffect, useMemo, useState } from 'react';
import { Drawer, Toolbar } from '@material-ui/core';
import MapGL, { MapEvent } from 'react-map-gl';
import { MapboxGeoJSONFeature } from 'mapbox-gl';

import { FeatureSidebar } from './FeatureSidebar';
import { MapTooltip } from './map/MapTooltip';
import { MapParams, useMapContent } from './map/use-map-content';
import { BackgroundControl } from './controls/BackgroundControl';
import { NetworkControl } from './controls/NetworkControl';
import { useLayerSelection } from './controls/use-layer-selection';
import { ViewName, views } from './config/views';
import { BackgroundName } from './config/types';
import { LayerName, layers } from './config/layers';
import { MapStyle } from './map/MapStyle';
import { HazardsControl } from './controls/HazardsControl';

const viewportLimits = {
  minZoom: 3,
  maxZoom: 16,
  maxPitch: 0,
};

const MAPBOX_KEY = 'pk.eyJ1IjoidG9tcnVzc2VsbCIsImEiOiJjaXZpMTFpdGkwMDQ1MnptcTh4ZzRzeXNsIn0.ZSvSOHSsWBQ44QNhA71M6Q';

interface MapViewProps {
  view: ViewName;
}

export const MapView: FC<MapViewProps> = ({ view }) => {
  const [viewport, setViewport] = useState({
    latitude: 18.14,
    longitude: -77.28,
    zoom: 8,
  });

  const [background, setBackground] = useState<BackgroundName>('light');

  const [hoveredFeatures, setHoveredFeatures] = useState<MapboxGeoJSONFeature[]>([]);
  const [hoverLngLat, setHoverLngLat] = useState<[number, number]>(null);
  const handleMapHover = useCallback((e: MapEvent) => {
    setHoveredFeatures(e.features ?? []);
    if (e.features?.length) {
      setHoverLngLat(e.lngLat);
    }
  }, []);

  const [selectedFeatures, setSelectedFeatures] = useState<MapboxGeoJSONFeature[]>([]);
  const handleMapClick = useCallback((e: MapEvent) => {
    setSelectedFeatures(e.features ?? []);
  }, []);

  const viewLayerNames = useMemo<LayerName[]>(() => views[view].layers as LayerName[], [view]);
  const layerDefinitions = useMemo(
    () => viewLayerNames.map((layerName) => ({ ...layers[layerName], key: layerName })),
    [viewLayerNames],
  );

  const { layerSelection, updateLayerSelection, selectSingleLayer } = useLayerSelection(viewLayerNames);

  // TODO: create separate mechanism for layer-dependent features
  const [fluvialReturnPeriod, setFluvialReturnPeriod] = useState(10);
  useEffect(() => {
    if (view === 'hazards') {
      selectSingleLayer(`flood_fluvial_${fluvialReturnPeriod}` as LayerName);
    }
  }, [view, fluvialReturnPeriod, selectSingleLayer]);

  const [coastalReturnPeriod, setCoastalReturnPeriod] = useState(1);
  useEffect(() => {
    if (view === 'hazards') {
      selectSingleLayer(`flood_coastal_${coastalReturnPeriod}` as LayerName);
    }
  }, [view, coastalReturnPeriod, selectSingleLayer]);

  const mapContentParams = useMemo<MapParams>(
    () => ({
      background,
      view,
      dataLayerSelection: layerSelection,
      highlightedFeature: selectedFeatures?.[0],
    }),
    [background, view, layerSelection, selectedFeatures],
  );
  const mapContent = useMapContent(mapContentParams);

  return (
    <>
      <Drawer variant="permanent">
        <Toolbar /> {/* Prevents app bar from concealing content*/}
        <div className="drawer-contents">
          {view === 'overview' && (
            <NetworkControl
              dataLayers={layerDefinitions}
              layerVisibility={layerSelection}
              onLayerVisChange={updateLayerSelection}
            />
          )}
          {view === 'hazards' && (
            <HazardsControl
              fluvialReturnPeriod={fluvialReturnPeriod}
              setFluvialReturnPeriod={setFluvialReturnPeriod}
              coastalReturnPeriod={coastalReturnPeriod}
              setCoastalReturnPeriod={setCoastalReturnPeriod}
            />
          )}
          <BackgroundControl background={background} onBackgroundChange={setBackground} />
        </div>
      </Drawer>
      <div className="map-height">
        <MapGL
          mapboxApiAccessToken={MAPBOX_KEY}
          width="100%"
          height="100%"
          {...viewport}
          onViewportChange={setViewport}
          {...viewportLimits}
          dragRotate={false}
          touchRotate={false}
          onHover={handleMapHover}
          onClick={handleMapClick}
          reuseMaps={true}
        >
          <MapStyle mapStyle={mapContent} />
          {hoveredFeatures.length !== 0 && <MapTooltip features={hoveredFeatures} tooltipLngLat={hoverLngLat} />}
        </MapGL>
        {selectedFeatures.length !== 0 && <FeatureSidebar feature={selectedFeatures[0]} />}
      </div>
    </>
  );
};
