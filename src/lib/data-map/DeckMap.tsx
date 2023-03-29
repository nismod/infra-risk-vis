import DeckGL, { DeckGLContextValue, DeckGLRef, DeckProps, MapView, MapViewState } from 'deck.gl/typed';
import { FC, Provider, createContext, useRef, useState } from 'react';

import { useTriggerMemo } from '../hooks/use-trigger-memo';
import { MapContextProviderWithLimits } from './MapContextProviderWithLimits';

interface DeckMapProps {
  initialViewState: any;
  layersFunction: ({ zoom }) => any[];
  dataLoadTrigger?: number;
  onHover: any;
  onClick?: any;
  layerRenderFilter: DeckProps['layerFilter'];
  pickingRadius?: number;
}

export const ViewStateContext = createContext<{
  viewState: MapViewState;
  setViewState: (viewState: MapViewState) => void;
}>(null);

export const DeckMap: FC<DeckMapProps> = ({
  initialViewState,
  layersFunction,
  dataLoadTrigger,
  onHover,
  onClick,
  layerRenderFilter,
  pickingRadius,
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
        views={[
          new MapView({
            repeat: true,
            controller: {
              scrollZoom: {
                smooth: true,
                speed: 0.2,
              },
              keyboard: false, //can't deactivate keyboard rotate only so deactivate all keyboard
              dragRotate: false,
              touchRotate: false,
            },
          }),
        ]}
        viewState={viewState}
        onViewStateChange={({ viewState }) => setViewState(viewState)}
        layers={layers}
        layerFilter={layerRenderFilter}
        onHover={(info, event) =>
          !event.srcEvent.defaultPrevented && // ignore pointer events from HUD: https://github.com/visgl/deck.gl/discussions/6252
          deckRef.current &&
          onHover(info, deckRef.current)
        }
        onClick={(info, event) =>
          !event.srcEvent.defaultPrevented && // ignore pointer events from HUD: https://github.com/visgl/deck.gl/discussions/6252
          deckRef.current &&
          onClick?.(info, deckRef.current)
        }
        pickingRadius={pickingRadius}
        ContextProvider={
          MapContextProviderWithLimits as unknown as Provider<DeckGLContextValue> /* unknown because TS doesn't like the cast */
        }
      >
        {/* make sure components like StaticMap are immediate children of DeckGL so that they 
            can be managed properly by Deck - see https://deck.gl/docs/api-reference/react/deckgl#jsx-layers
        */}
        {children}
      </DeckGL>
    </ViewStateContext.Provider>
  );
};
