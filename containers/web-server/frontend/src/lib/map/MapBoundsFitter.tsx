import { easeCubic } from 'd3-ease';
import { FC, useContext } from 'react';
import { FlyToInterpolator, WebMercatorViewport } from 'react-map-gl';

import { BoundingBox, appToDeckBoundingBox } from '@/lib/bounding-box';
import { ViewStateContext } from '@/lib/data-map/DeckMap';
import { useChangeEffect } from '@/lib/hooks/use-change-effect';

interface MapBoundsFitterProps {
  boundingBox: BoundingBox;
}

export const MapBoundsFitter: FC<MapBoundsFitterProps> = ({ boundingBox }) => {
  const { viewState, setViewState } = useContext(ViewStateContext);

  useChangeEffect(
    () => {
      if (boundingBox != null) {
        const deckBbox = appToDeckBoundingBox(boundingBox);
        const viewport = new WebMercatorViewport({ width: 800, height: 600 });
        const { latitude, longitude, zoom } = viewport.fitBounds(deckBbox, { padding: 20 });

        setViewState({
          ...viewState,
          latitude,
          longitude,
          zoom,
          transitionDuration: 1500,
          transitionInterpolator: new FlyToInterpolator() as any,
          transitionEasing: easeCubic,
        });
      }
    },
    [boundingBox, viewState, setViewState],
    [boundingBox],
  );

  return null;
};
