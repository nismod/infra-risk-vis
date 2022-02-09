import DeckGL, { DeckProps } from 'deck.gl';
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
  layersFunction: ({ zoom }) => any[];
  backgroundStyle: object;
  onHover: any;
  onClick?: any;
  layerRenderFilter: DeckProps['layerFilter'];
  pickingRadius?: number;
  targetBounds: BoundingBox;
}

export const MapViewport: FC<MapViewportProps> = ({
  initialViewState,
  layersFunction,
  backgroundStyle,
  onHover,
  onClick,
  layerRenderFilter,
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
      layerFilter={layerRenderFilter}
      onHover={(info) => deckRef.current && onHover(info, deckRef.current)}
      onClick={(info) => deckRef.current && onClick?.(info, deckRef.current)}
      pickingRadius={pickingRadius}
      ContextProvider={MapContext.Provider}
    >
      <StaticMap mapStyle={backgroundStyle} attributionControl={false} />
      <MapBoundsFitter boundingBox={targetBounds} viewState={viewState} onViewState={setViewState} />

      {children}
    </DeckGL>
  );
};
