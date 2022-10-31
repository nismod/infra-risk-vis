import _ from 'lodash';
import { FC, ReactNode, useCallback, useEffect, useMemo } from 'react';
import { StaticMap } from 'react-map-gl';

import { usePrevious } from '@/lib/hooks/use-previous';
import { useTrackingRef } from '@/lib/hooks/use-tracking-ref';
import { useTrigger } from '@/lib/hooks/use-trigger';

import { DeckMap } from './DeckMap';
import { useInteractions } from './interactions/use-interactions';
import { ViewLayer, ViewLayerParams } from './view-layers';

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

  const dataLoaders = useMemo(
    () => viewLayers.map((vl) => vl.dataAccessFn?.(vl.styleParams?.colorMap?.fieldSpec)?.dataLoader).filter(Boolean),
    [viewLayers],
  );

  const [dataLoadTrigger, triggerDataUpdate] = useTrigger();

  const doTrigger = useCallback(() => {
    triggerDataUpdate();
  }, [triggerDataUpdate]);

  const previousLoaders = usePrevious(dataLoaders);

  useEffect(() => {
    // destroy removed data loaders to free up memory
    const removedLoaders = _.difference(previousLoaders ?? [], dataLoaders);
    removedLoaders.forEach((dl) => dl.destroy());

    // subscribe to new data loaders to get notified when data is loaded
    const addedLoaders = _.difference(dataLoaders, previousLoaders ?? []);
    addedLoaders.forEach((dl) => dl.subscribe(doTrigger));

    // if there was a change in data loaders, trigger an update to the data map
    if (addedLoaders.length > 0 || removedLoaders.length > 0) {
      doTrigger();
    }
  }, [dataLoaders, previousLoaders, doTrigger]);

  /* store current value of dataLoaders so that we can clean up data on component unmount
   * this is necessary because we don't want to keep the data loaders around after the component is unmounted
   */
  const currentLoadersRef = useTrackingRef(dataLoaders);
  useEffect(() => {
    return () => {
      // eslint-disable-next-line react-hooks/exhaustive-deps
      currentLoadersRef.current?.forEach((dl) => dl.destroy());
    };
  }, [currentLoadersRef]);

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
