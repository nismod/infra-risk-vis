import DeckGL from 'deck.gl';
import { easeCubic } from 'd3-ease';
import { FC, useEffect, useMemo, useRef, useState } from 'react';
import {
  AttributionControl,
  MapContext,
  MapContextProps,
  NavigationControl,
  ScaleControl,
  StaticMap,
  FlyToInterpolator,
} from 'react-map-gl';
import _ from 'lodash';

import { backgroundConfig, BackgroundName } from '../config/backgrounds';

function visible(isVisible: boolean): 'visible' | 'none' {
  return isVisible ? 'visible' : 'none';
}

function makeMapboxConfig(background: BackgroundName) {
  return {
    version: 8,
    sources: _.mapValues(backgroundConfig, (b) => b.source),
    layers: Object.values(backgroundConfig).map((b) =>
      _.merge(b.layer, { layout: { visibility: visible(background === b.id) } }),
    ),
  };
}

interface ViewTarget {
  latitude: number;
  longitude: number;
  zoom: number;
}

interface MapViewportProps {
  layersFunction: any;
  background: BackgroundName;
  onHover: any;
  onClick: any;
  pickingRadius: number;
  viewTarget: ViewTarget;
}

export const MapViewport: FC<MapViewportProps> = ({
  layersFunction,
  background,
  onHover,
  onClick,
  pickingRadius,
  viewTarget,
  children,
}) => {
  const [viewport, setViewport] = useState<any>({
    latitude: 18.14,
    longitude: -77.28,
    zoom: 8,
    minZoom: 3,
    maxZoom: 16,
    maxPitch: 0,
  });

  useEffect(() => {
    if (viewTarget != null) {
      setViewport({
        ...viewport,
        ...viewTarget,
        transitionDuration: 1500,
        transitionInterpolator: new FlyToInterpolator(),
        transitionEasing: easeCubic,
      });
    }
  }, [viewTarget]);

  const deckRef = useRef<DeckGL<MapContextProps>>();

  const zoom = viewport.zoom;

  const backgroundStyle = useMemo(() => makeMapboxConfig(background), [background]);
  const layers = useMemo(() => layersFunction({ zoom }), [layersFunction, zoom]);

  return (
    <DeckGL<MapContextProps>
      ref={deckRef}
      style={{
        overflow: 'hidden',
      }}
      getCursor={() => 'default'}
      controller={true}
      viewState={viewport}
      onViewStateChange={({ viewState }) => setViewport(viewState)}
      layers={layers}
      layerFilter={({ layer, renderPass }) => {
        if (renderPass === 'picking:hover') {
          // don't render raster layers on hover picking pass (but render them for manual picking)
          if (layer.id.match(/^(coastal|fluvial|surface|cyclone)/)) return false;
        }
        return true;
      }}
      pickingRadius={pickingRadius}
      onHover={(info) => deckRef.current && onHover(info, deckRef.current)}
      onClick={(info) => deckRef.current && onClick(info, deckRef.current)}
      ContextProvider={MapContext.Provider}
    >
      <StaticMap mapStyle={backgroundStyle} attributionControl={false} />
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
      {children}
    </DeckGL>
  );
};
