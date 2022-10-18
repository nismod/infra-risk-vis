import { Box } from '@mui/material';
import DeckGL, { DeckProps } from 'deck.gl';
import { FC, ReactNode, createContext, useMemo, useRef, useState } from 'react';
import { MapContext, MapContextProps } from 'react-map-gl';

interface DeckMapProps {
  initialViewState: any;
  layersFunction: ({ zoom }) => any[];
  dataLoadTrigger?: number;
  onHover: any;
  onClick?: any;
  layerRenderFilter: DeckProps['layerFilter'];
  pickingRadius?: number;
  uiOverlays: ReactNode;
}

export const ViewStateContext = createContext<{
  viewState: any;
  setViewState: (viewState: any) => void;
}>(null);

export const DeckMap: FC<DeckMapProps> = ({
  initialViewState,
  layersFunction,
  dataLoadTrigger,
  onHover,
  onClick,
  layerRenderFilter,
  pickingRadius,
  uiOverlays,
  children,
}) => {
  const [viewState, setViewState] = useState<any>(initialViewState);

  const deckRef = useRef<DeckGL<MapContextProps>>();

  const zoom = viewState.zoom;

  const layers = useMemo(() => layersFunction({ zoom }), [layersFunction, zoom, dataLoadTrigger]);

  return (
    <ViewStateContext.Provider value={{ viewState, setViewState }}>
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
        {children}
      </DeckGL>
      {uiOverlays && (
        <div
          style={{
            pointerEvents: 'none',
            position: 'absolute',
            top: 0,
            left: 0,
            width: '100%',
            height: '100%',
            zIndex: 1,
          }}
        >
          <Box sx={{ pointerEvents: 'auto' }}>{uiOverlays}</Box>
        </div>
      )}
    </ViewStateContext.Provider>
  );
};
