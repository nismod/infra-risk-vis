import React, { FC, useCallback, useMemo, useState } from 'react';
import { Drawer, Toolbar } from '@material-ui/core';
import { StaticMap } from 'react-map-gl';
import { MapboxGeoJSONFeature } from 'mapbox-gl';
import { readPixelsToArray } from '@luma.gl/core';

import DeckGL from 'deck.gl';

import { FeatureSidebar } from './FeatureSidebar';
import { MapTooltip } from './map/tooltip/MapTooltip';
import { FeaturesTooltipContent } from './map/tooltip/FeaturesTooltipContent';
import { MapParams, useDeckLayers, useMapboxStyle } from './map/use-map-content';
import { BackgroundControl } from './controls/BackgroundControl';
import { NetworkControl } from './controls/NetworkControl';
import { useLayerSelection } from './controls/use-layer-selection';
import { ViewName, views } from './config/views';
import { BackgroundName } from './config/types';
import { LayerName, layers } from './config/layers';
import { HazardsControl } from './controls/HazardsControl';

const MAPBOX_KEY = 'pk.eyJ1IjoidG9tcnVzc2VsbCIsImEiOiJjaXZpMTFpdGkwMDQ1MnptcTh4ZzRzeXNsIn0.ZSvSOHSsWBQ44QNhA71M6Q';

interface MapViewProps {
  view: ViewName;
}

export const MapView: FC<MapViewProps> = ({ view }) => {
  const [viewport, setViewport] = useState({
    latitude: 18.14,
    longitude: -77.28,
    zoom: 8,
    minZoom: 3,
    maxZoom: 16,
    maxPitch: 0,
  });

  const [background, setBackground] = useState<BackgroundName>('light');

  const [hoverColor, setHoverColor] = useState<[number, number, number, number]>(null);
  const [hoveredFeatures, setHoveredFeatures] = useState<MapboxGeoJSONFeature[]>([]);
  const [hoverXY, setHoverXY] = useState<[number, number]>(null);

  const [selectedFeatures, setSelectedFeatures] = useState<MapboxGeoJSONFeature[]>([]);

  const viewLayerNames = useMemo<LayerName[]>(() => views[view].layers as LayerName[], [view]);
  const layerDefinitions = useMemo(
    () => viewLayerNames.map((layerName) => ({ ...layers[layerName], key: layerName })),
    [viewLayerNames],
  );

  const { layerSelection, updateLayerSelection } = useLayerSelection(viewLayerNames);

  const mapContentParams = useMemo<MapParams>(
    () => ({
      background,
      view,
      dataLayerSelection: layerSelection,
      highlightedFeature: selectedFeatures?.[0],
    }),
    [background, view, layerSelection, selectedFeatures],
  );

  const onLayerHover = useCallback((info: any) => {
    const { bitmap, object, sourceLayer, x, y } = info;

    if (bitmap || object) {
      setHoverXY([x, y]);
    } else {
      setHoverXY(null);
    }
    if (bitmap) {
      const pixelColor = readPixelsToArray(sourceLayer.props.image, {
        sourceX: bitmap.pixel[0],
        sourceY: bitmap.pixel[1],
        sourceWidth: 1,
        sourceHeight: 1,
        sourceType: undefined,
      });
      if (pixelColor[3]) {
        setHoverColor(pixelColor);
      } else {
        setHoverColor(null);
      }
      setHoveredFeatures(null);
    } else if (object?.properties) {
      setHoverXY([x, y]);
      setHoveredFeatures([object]);
    } else {
      setHoverColor(null);
      setHoveredFeatures(null);
    }
  }, []);

  const onLayerClick = useCallback((info: any) => {
    if (!info.bitmap && info.object) {
      setSelectedFeatures([info.object]);
    } else {
      setSelectedFeatures([]);
    }
  }, []);

  const mapboxStyle = useMapboxStyle(mapContentParams);
  const deckLayers = useDeckLayers(mapContentParams, onLayerHover);

  return (
    <>
      <Drawer variant="permanent">
        <Toolbar /> {/* Prevents app bar from concealing content*/}
        <div className="drawer-contents">
          {view === 'overview' && (
            <>
              <NetworkControl
                dataLayers={layerDefinitions}
                layerVisibility={layerSelection}
                onLayerVisChange={updateLayerSelection}
              />
              <HazardsControl layerVisibility={layerSelection} onLayerVisibilityUpdate={updateLayerSelection} />
            </>
          )}
          <BackgroundControl background={background} onBackgroundChange={setBackground} />
        </div>
      </Drawer>
      <div className="map-height">
        <DeckGL
          controller={true}
          viewState={viewport}
          onViewStateChange={({ viewState }) => setViewport(viewState)}
          layers={deckLayers}
          pickingRadius={10}
          onHover={onLayerHover}
          onClick={onLayerClick}
        >
          <StaticMap mapStyle={mapboxStyle} mapboxApiAccessToken={MAPBOX_KEY} />
          <MapTooltip tooltipXY={hoverXY}>
            {hoverColor && (
              <div style={{ backgroundColor: `rgb(${hoverColor[0]},${hoverColor[1]},${hoverColor[2]})` }}>Flooding</div>
            )}
            {hoveredFeatures && <FeaturesTooltipContent features={hoveredFeatures} />}
          </MapTooltip>
        </DeckGL>
        {selectedFeatures.length !== 0 && <FeatureSidebar feature={selectedFeatures[0]} />}
      </div>
    </>
  );
};
