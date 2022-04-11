import { FC, ReactNode, useCallback, useEffect } from 'react';
import { StaticMap } from 'react-map-gl';

import { ViewLayer, ViewLayerParams } from './view-layers';

import { DeckMap } from './DeckMap';
import { useInteractions } from './interactions/use-interactions';
import { useTrigger } from 'lib/hooks/use-trigger';

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

  const [dataLoadTrigger, triggerDataUpdate] = useTrigger();

  useEffect(() => {
    const dataLoaders = viewLayers.map((vl) => vl.dataAccessFn?.(viewLayersParams[vl.id])?.dataLoader).filter(Boolean);

    dataLoaders.forEach((dl) => dl.subscribe(triggerDataUpdate));

    return () => {
      dataLoaders.forEach((dl) => dl.unsubscribe(triggerDataUpdate));
    };
  }, [viewLayers, viewLayersParams, triggerDataUpdate]);

  const deckLayersFunction = useCallback(
    ({ zoom }: { zoom: number }) =>
      viewLayers.map((viewLayer) => makeDeckLayers(viewLayer, viewLayersParams[viewLayer.id], zoom)),
    [viewLayers, viewLayersParams],
  );

  return (
    <DeckMap
      initialViewState={initialViewState}
      layersFunction={deckLayersFunction}
      dataLoadTrigger={dataLoadTrigger}
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
