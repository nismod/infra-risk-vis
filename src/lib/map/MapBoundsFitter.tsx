import { useContext } from 'react';
import { MapContext, FlyToInterpolator } from 'react-map-gl';

import { easeCubic } from 'd3-ease';

import { useChangeEffect } from 'lib/hooks/use-change-effect';

type DeckBoundingBox = [[number, number], [number, number]];

export const MapBoundsFitter = ({ boundingBox, viewState, onViewState }) => {
  const { viewport } = useContext(MapContext);

  useChangeEffect(
    () => {
      if (boundingBox != null) {
        const deckBbox: DeckBoundingBox = [
          [boundingBox.minX, boundingBox.minY],
          [boundingBox.maxX, boundingBox.maxY],
        ];
        const { latitude, longitude, zoom } = viewport.fitBounds(deckBbox, { padding: 20 });
        onViewState({
          ...viewState,
          latitude,
          longitude,
          zoom,
          transitionDuration: 1500,
          transitionInterpolator: new FlyToInterpolator(),
          transitionEasing: easeCubic,
        });
      }
    },
    [boundingBox, viewState, onViewState, viewport],
    [boundingBox],
  );

  return null;
};
