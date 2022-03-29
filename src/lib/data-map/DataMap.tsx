import { FC, ReactNode, useCallback } from 'react';
import { StaticMap } from 'react-map-gl';

import { ViewLayer, ViewLayerParams } from './view-layers';

import { DeckMap } from './DeckMap';
import { useInteractions } from './interactions/use-interactions';

export interface DataMapProps {
  initialViewState: any;
  viewLayers: ViewLayer[];
  viewLayersParams: Record<string, ViewLayerParams>;
  interactionGroups: any;
  backgroundStyle: any;
  uiOverlays?: ReactNode;
}

// set a convention where the view layer id is either the first part of the deck id before the @ sign, or it's the whole id
function lookupViewForDeck(deckLayerId: string) {
  return deckLayerId.split('@')[0];
}

export const DataMap: FC<DataMapProps> = ({
  initialViewState,
  viewLayers,
  viewLayersParams,
  interactionGroups,
  backgroundStyle,
  uiOverlays,
  children,
}) => {
  const { onHover, onClick, layerFilter, pickingRadius } = useInteractions(
    viewLayers,
    lookupViewForDeck,
    interactionGroups,
  );

  const deckLayersFunction = useCallback(
    ({ zoom }: { zoom: number }) =>
      viewLayers.map((viewLayer) => makeDeckLayers(viewLayer, viewLayersParams[viewLayer.id], zoom)),
    [viewLayers, viewLayersParams],
  );

  return (
    <DeckMap
      initialViewState={initialViewState}
      layersFunction={deckLayersFunction}
      onHover={onHover}
      onClick={onClick}
      layerRenderFilter={layerFilter}
      pickingRadius={pickingRadius}
      uiOverlays={uiOverlays}
    >
      <StaticMap mapStyle={backgroundStyle} attributionControl={false} />
      {children}
    </DeckMap>
  );
};

function makeDeckLayers(viewLayer: ViewLayer, viewLayerParams: ViewLayerParams, zoom: number) {
  return viewLayer.fn({
    deckProps: { id: viewLayer.id, pickable: !!viewLayer.interactionGroup },
    zoom,
    ...viewLayerParams,
  });
}
