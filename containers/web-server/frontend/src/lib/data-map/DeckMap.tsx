import { Box } from '@mui/material';
import DeckGL, { DeckGLContextValue, DeckGLRef, DeckProps } from 'deck.gl/typed';
import { FC, Provider, ReactNode, createContext, useRef, useState } from 'react';
import { MapContext } from 'react-map-gl';

import { useTriggerMemo } from '../hooks/use-trigger-memo';

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

  const deckRef = useRef<DeckGLRef>();

  const zoom = viewState.zoom;

  const layers = useTriggerMemo(() => layersFunction({ zoom }), [layersFunction, zoom], dataLoadTrigger);

  return (
    <ViewStateContext.Provider value={{ viewState, setViewState }}>
      <DeckGL
        ref={deckRef}
        style={{
          overflow: 'hidden',
        }}
        getCursor={() => 'default'}
        controller={{
          keyboard: false, //can't deactivate keyboard rotate only so deactivate all keyboard
          dragRotate: false,
          touchRotate: false,
        }}
        viewState={viewState}
        onViewStateChange={({ viewState }) => setViewState(viewState)}
        layers={layers}
        layerFilter={layerRenderFilter}
        onHover={(info) => deckRef.current && onHover(info, deckRef.current)}
        onClick={(info) => deckRef.current && onClick?.(info, deckRef.current)}
        pickingRadius={pickingRadius}
        ContextProvider={
          MapContext.Provider as unknown as Provider<DeckGLContextValue> /* unknown because TS doesn't like the cast */
        }
      >
        {/* make sure components like StaticMap are immediate children of DeckGL so that they 
            can be managed properly by Deck - see https://deck.gl/docs/api-reference/react/deckgl#jsx-layers
        */}
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
