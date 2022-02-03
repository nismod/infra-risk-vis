import DeckGL from 'deck.gl';
import { FC, useMemo, useRef, useState } from 'react';
import { MapContext, MapContextProps, StaticMap } from 'react-map-gl';

import { MapBoundsFitter } from './MapBoundsFitter';

export interface BoundingBox {
  minX: number;
  minY: number;
  maxX: number;
  maxY: number;
}

interface MapViewportProps {
  initialViewState: any;
  layersFunction: any;
  backgroundStyle: object;
  onHover: any;
  onClick: any;
  pickingRadius: number;
  targetBounds: BoundingBox;
  attribution?: string;
}

export const MapViewport: FC<MapViewportProps> = ({
  initialViewState,
  layersFunction,
  backgroundStyle,
  onHover,
  onClick,
  pickingRadius,
  targetBounds,
  children,
}) => {
  const [viewState, setViewState] = useState<any>(initialViewState);

  const deckRef = useRef<DeckGL<MapContextProps>>();

  const zoom = viewState.zoom;

  const layers = useMemo(() => layersFunction({ zoom }), [layersFunction, zoom]);

  return (
    <DeckGL<MapContextProps>
      ref={deckRef}
      style={{
        overflow: 'hidden',
      }}
      getCursor={() => 'default'}
      controller={true}
      viewState={viewState}
      onViewStateChange={({ viewState }) => setViewState(viewState)}
      layers={layers}
      layerFilter={({ layer, renderPass }) => {
        if (renderPass === 'picking:hover') {
          // don't render raster and region layers on hover picking pass (but render them for manual picking)
          if (layer.id.match(/^(coastal|fluvial|surface|cyclone|boundaries)/)) return false;
        }
        return true;
      }}
      pickingRadius={pickingRadius}
      onHover={(info) => deckRef.current && onHover(info, deckRef.current)}
      onClick={(info) => deckRef.current && onClick(info, deckRef.current)}
      ContextProvider={MapContext.Provider}
    >
      <StaticMap mapStyle={backgroundStyle} attributionControl={false} />
      <MapBoundsFitter boundingBox={targetBounds} viewState={viewState} onViewState={setViewState} />

      {children}
    </DeckGL>
  );
};
